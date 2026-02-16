"""
Performance Optimization Module for EPyR Tools
===============================================

This module provides performance optimizations for handling large EPR datasets:
- Memory-efficient data loading
- Chunked processing for large files
- Caching system for frequently accessed data
- Memory usage monitoring and optimization

Usage:
    from epyr.performance import OptimizedLoader, DataCache

    # Use optimized loader for large files
    loader = OptimizedLoader(chunk_size_mb=10, cache_enabled=True)
    x, y, params = loader.load_epr_file(file_path)

    # Use caching for repeated access
    cache = DataCache(max_size_mb=100)
    cached_data = cache.get_or_load(file_path, load_function)
"""

import gc
import hashlib
import os
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from .config import config
from .logging_config import get_logger

logger = get_logger(__name__)


class MemoryMonitor:
    """Monitor and optimize memory usage during data processing."""

    @staticmethod
    def get_memory_info() -> Dict[str, float]:
        """Get current memory usage information.

        Returns:
            Dict with memory info in MB: {rss, vms, percent}
        """
        try:
            import psutil

            process = psutil.Process()
            memory_info = process.memory_info()
            memory_percent = process.memory_percent()

            return {
                "rss_mb": memory_info.rss / (1024 * 1024),  # Resident Set Size
                "vms_mb": memory_info.vms / (1024 * 1024),  # Virtual Memory Size
                "percent": memory_percent,
            }
        except ImportError:
            logger.warning("psutil not available for memory monitoring")
            return {"rss_mb": 0, "vms_mb": 0, "percent": 0}

    @staticmethod
    def check_memory_limit() -> bool:
        """Check if memory usage is approaching configured limit.

        Returns:
            True if memory usage is acceptable, False if limit exceeded
        """
        memory_limit_mb = config.get("performance.memory_limit_mb", 500)
        memory_info = MemoryMonitor.get_memory_info()

        if memory_info["rss_mb"] > memory_limit_mb:
            logger.warning(
                f"Memory usage ({memory_info['rss_mb']:.1f} MB) exceeds "
                f"configured limit ({memory_limit_mb} MB)"
            )
            return False
        return True

    @staticmethod
    def optimize_memory():
        """Perform memory optimization steps."""
        logger.debug("Running garbage collection")
        collected = gc.collect()
        if collected > 0:
            logger.debug(f"Garbage collector freed {collected} objects")


class DataCache:
    """LRU cache for frequently accessed EPR data files."""

    def __init__(self, max_size_mb: Optional[int] = None):
        """Initialize data cache.

        Args:
            max_size_mb: Maximum cache size in MB. Uses config default if None.
        """
        self.max_size_mb = max_size_mb or config.get("performance.cache_size_mb", 100)
        self.cache: Dict[str, Dict[str, Any]] = {}
        self.access_order: List[str] = []

    def _get_file_hash(self, file_path: Path) -> str:
        """Generate cache key from file path and modification time."""
        stat = file_path.stat()
        hash_str = f"{file_path}_{stat.st_mtime}_{stat.st_size}"
        return hashlib.md5(hash_str.encode()).hexdigest()

    def _get_data_size_mb(self, data: Tuple) -> float:
        """Estimate memory size of cached data in MB."""
        total_bytes = 0
        for item in data:
            if isinstance(item, np.ndarray):
                total_bytes += item.nbytes
            elif isinstance(item, (list, tuple)):
                for sub_item in item:
                    if isinstance(sub_item, np.ndarray):
                        total_bytes += sub_item.nbytes
            elif isinstance(item, dict):
                # Rough estimate for parameter dictionaries
                total_bytes += len(str(item).encode("utf-8"))
        return total_bytes / (1024 * 1024)

    def _evict_if_needed(self, new_data_size_mb: float):
        """Evict least recently used items if cache would exceed limit."""
        current_size = sum(entry["size_mb"] for entry in self.cache.values())

        while current_size + new_data_size_mb > self.max_size_mb and self.access_order:
            oldest_key = self.access_order.pop(0)
            if oldest_key in self.cache:
                evicted_size = self.cache[oldest_key]["size_mb"]
                del self.cache[oldest_key]
                current_size -= evicted_size
                logger.debug(f"Evicted cached data ({evicted_size:.1f} MB)")

    def get(self, file_path: Path) -> Optional[Tuple]:
        """Get cached data for file if available and still valid."""
        if not config.get("performance.cache_enabled", True):
            return None

        cache_key = self._get_file_hash(file_path)

        if cache_key in self.cache:
            # Move to end of access order (most recently used)
            self.access_order.remove(cache_key)
            self.access_order.append(cache_key)
            logger.debug(f"Cache hit for {file_path.name}")
            return self.cache[cache_key]["data"]

        return None

    def put(self, file_path: Path, data: Tuple):
        """Cache data for file."""
        if not config.get("performance.cache_enabled", True):
            return

        cache_key = self._get_file_hash(file_path)
        data_size_mb = self._get_data_size_mb(data)

        # Don't cache if single item exceeds cache limit
        if data_size_mb > self.max_size_mb:
            logger.debug(f"Data too large to cache ({data_size_mb:.1f} MB)")
            return

        self._evict_if_needed(data_size_mb)

        self.cache[cache_key] = {
            "data": data,
            "size_mb": data_size_mb,
        }
        self.access_order.append(cache_key)
        logger.debug(f"Cached data for {file_path.name} ({data_size_mb:.1f} MB)")

    def clear(self):
        """Clear all cached data."""
        self.cache.clear()
        self.access_order.clear()
        logger.debug("Cache cleared")

    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        total_size_mb = sum(entry["size_mb"] for entry in self.cache.values())
        return {
            "entries": len(self.cache),
            "total_size_mb": total_size_mb,
            "max_size_mb": self.max_size_mb,
            "usage_percent": (
                (total_size_mb / self.max_size_mb) * 100 if self.max_size_mb > 0 else 0
            ),
        }


class OptimizedLoader:
    """Optimized data loader for large EPR datasets."""

    def __init__(self, chunk_size_mb: Optional[int] = None, cache_enabled: bool = True):
        """Initialize optimized loader.

        Args:
            chunk_size_mb: Chunk size for processing large files
            cache_enabled: Whether to use caching
        """
        self.chunk_size_mb = chunk_size_mb or config.get(
            "performance.chunk_size_mb", 10
        )
        self.cache = DataCache() if cache_enabled else None

    def load_epr_file(self, file_path: Union[str, Path]) -> Tuple:
        """Load EPR file with optimization for large datasets.

        Args:
            file_path: Path to EPR file

        Returns:
            Tuple of (x_data, y_data, parameters, file_path_str)
        """
        file_path = Path(file_path)

        # Check cache first
        if self.cache:
            cached_data = self.cache.get(file_path)
            if cached_data is not None:
                return cached_data

        # Check memory before loading
        if not MemoryMonitor.check_memory_limit():
            MemoryMonitor.optimize_memory()

        # Load data using existing eprload function
        from .eprload import eprload

        logger.debug(f"Loading {file_path} with optimization")

        try:
            # Load with plotting disabled for performance
            data = eprload(str(file_path), plot_if_possible=False)

            # Cache the result
            if self.cache:
                self.cache.put(file_path, data)

            return data

        except Exception as e:
            logger.error(f"Optimized loading failed for {file_path}: {e}")
            raise

    def load_chunked_data(
        self, file_path: Union[str, Path], chunk_processor: Callable
    ) -> Any:
        """Load and process large data files in chunks.

        Args:
            file_path: Path to data file
            chunk_processor: Function to process each chunk

        Returns:
            Processed result from chunk_processor
        """
        file_path = Path(file_path)
        file_size_mb = file_path.stat().st_size / (1024 * 1024)

        if file_size_mb <= self.chunk_size_mb:
            # File is small enough to load entirely
            return chunk_processor(file_path, is_chunked=False)

        logger.info(f"Processing large file ({file_size_mb:.1f} MB) in chunks")

        # For very large files, we would implement chunked reading
        # This is a framework for future implementation
        chunk_size_bytes = self.chunk_size_mb * 1024 * 1024

        # This is a placeholder - actual chunked implementation would depend
        # on the specific file format and processing requirements
        warnings.warn(
            f"Chunked processing not yet implemented for {file_path.suffix} files. "
            f"Loading entire file ({file_size_mb:.1f} MB)."
        )

        return chunk_processor(file_path, is_chunked=False)


def optimize_numpy_operations():
    """Configure NumPy for optimal performance."""
    # Set optimal number of threads for NumPy operations
    if config.get("performance.parallel_processing", True):
        try:
            import mkl

            # Use half of available cores to avoid oversubscription
            n_cores = max(1, os.cpu_count() // 2)
            mkl.set_num_threads(n_cores)
            logger.debug(f"Set MKL threads to {n_cores}")
        except ImportError:
            pass

        # Set OpenMP threads for NumPy
        os.environ["OMP_NUM_THREADS"] = str(max(1, os.cpu_count() // 2))


def get_performance_info() -> Dict[str, Any]:
    """Get comprehensive performance information.

    Returns:
        Dict with system performance metrics
    """
    memory_info = MemoryMonitor.get_memory_info()

    info = {
        "memory": memory_info,
        "cpu_count": os.cpu_count(),
        "config": {
            "cache_enabled": config.get("performance.cache_enabled"),
            "cache_size_mb": config.get("performance.cache_size_mb"),
            "chunk_size_mb": config.get("performance.chunk_size_mb"),
            "memory_limit_mb": config.get("performance.memory_limit_mb"),
            "parallel_processing": config.get("performance.parallel_processing"),
        },
    }

    return info


# Initialize performance optimizations on module import
optimize_numpy_operations()

# Global cache instance for convenience
_global_cache = DataCache()


def get_global_cache() -> DataCache:
    """Get the global data cache instance."""
    return _global_cache
