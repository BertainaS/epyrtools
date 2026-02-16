"""
Tests for EPyR Tools performance optimization module.
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

from epyr.performance import (
    DataCache,
    MemoryMonitor,
    OptimizedLoader,
    get_performance_info,
    optimize_numpy_operations,
)


class TestMemoryMonitor:
    """Test MemoryMonitor functionality."""

    @patch("epyr.performance.psutil")
    def test_get_memory_info_with_psutil(self, mock_psutil):
        """Test memory info retrieval with psutil available."""
        # Mock psutil
        mock_process = Mock()
        mock_process.memory_info.return_value = Mock(
            rss=100 * 1024 * 1024, vms=200 * 1024 * 1024  # 100 MB  # 200 MB
        )
        mock_process.memory_percent.return_value = 5.0
        mock_psutil.Process.return_value = mock_process

        memory_info = MemoryMonitor.get_memory_info()

        assert memory_info["rss_mb"] == 100.0
        assert memory_info["vms_mb"] == 200.0
        assert memory_info["percent"] == 5.0

    def test_get_memory_info_without_psutil(self):
        """Test memory info retrieval without psutil."""
        with patch("epyr.performance.psutil", None):
            with patch.dict("sys.modules", {"psutil": None}):
                memory_info = MemoryMonitor.get_memory_info()

                # Should return zero values when psutil not available
                assert memory_info["rss_mb"] == 0
                assert memory_info["vms_mb"] == 0
                assert memory_info["percent"] == 0

    @patch("epyr.performance.MemoryMonitor.get_memory_info")
    def test_check_memory_limit_ok(self, mock_get_memory):
        """Test memory limit check when within limits."""
        mock_get_memory.return_value = {"rss_mb": 100}  # 100 MB

        # Default limit is 500 MB
        result = MemoryMonitor.check_memory_limit()

        assert result is True

    @patch("epyr.performance.MemoryMonitor.get_memory_info")
    def test_check_memory_limit_exceeded(self, mock_get_memory):
        """Test memory limit check when limit exceeded."""
        mock_get_memory.return_value = {"rss_mb": 600}  # 600 MB

        # Default limit is 500 MB
        result = MemoryMonitor.check_memory_limit()

        assert result is False

    @patch("epyr.performance.gc.collect")
    def test_optimize_memory(self, mock_gc_collect):
        """Test memory optimization."""
        mock_gc_collect.return_value = 10  # 10 objects collected

        MemoryMonitor.optimize_memory()

        mock_gc_collect.assert_called_once()


class TestDataCache:
    """Test DataCache functionality."""

    def test_cache_initialization(self):
        """Test cache initialization."""
        cache = DataCache(max_size_mb=50)

        assert cache.max_size_mb == 50
        assert len(cache.cache) == 0
        assert len(cache.access_order) == 0

    def test_cache_initialization_from_config(self):
        """Test cache initialization using config values."""
        from epyr.config import config

        original_size = config.get("performance.cache_size_mb")
        try:
            config.set("performance.cache_size_mb", 75)
            cache = DataCache()
            assert cache.max_size_mb == 75
        finally:
            config.set("performance.cache_size_mb", original_size)

    def test_get_file_hash(self):
        """Test file hash generation."""
        cache = DataCache()

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test content")

        try:
            hash1 = cache._get_file_hash(test_file)
            hash2 = cache._get_file_hash(test_file)

            # Same file should generate same hash
            assert hash1 == hash2
            assert isinstance(hash1, str)
            assert len(hash1) == 32  # MD5 hash
        finally:
            test_file.unlink()

    def test_get_data_size_mb(self):
        """Test data size estimation."""
        cache = DataCache()

        # Test with numpy arrays
        x_data = np.array([1, 2, 3, 4, 5], dtype=np.float64)
        y_data = np.array([1, 4, 9, 16, 25], dtype=np.float64)
        metadata = {"key": "value"}

        data_tuple = (x_data, y_data, metadata)
        size_mb = cache._get_data_size_mb(data_tuple)

        assert size_mb > 0
        assert isinstance(size_mb, float)

        # Should be approximately the size of the arrays
        expected_size = (x_data.nbytes + y_data.nbytes) / (1024 * 1024)
        assert abs(size_mb - expected_size) < 0.001

    def test_cache_put_get(self):
        """Test basic cache put and get operations."""
        cache = DataCache(max_size_mb=10)

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test content")

        try:
            # Create test data
            test_data = (np.array([1, 2, 3]), np.array([4, 5, 6]), {"param": "value"})

            # Put data in cache
            cache.put(test_file, test_data)

            # Get data from cache
            retrieved_data = cache.get(test_file)

            assert retrieved_data is not None
            np.testing.assert_array_equal(retrieved_data[0], test_data[0])
            np.testing.assert_array_equal(retrieved_data[1], test_data[1])
            assert retrieved_data[2] == test_data[2]

        finally:
            test_file.unlink()

    def test_cache_miss(self):
        """Test cache miss for non-existent file."""
        cache = DataCache()

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)

        try:
            # File exists but no cached data
            result = cache.get(test_file)
            assert result is None
        finally:
            test_file.unlink()

    def test_cache_disabled(self):
        """Test cache behavior when disabled via config."""
        from epyr.config import config

        original_enabled = config.get("performance.cache_enabled")
        try:
            config.set("performance.cache_enabled", False)

            cache = DataCache()

            with tempfile.NamedTemporaryFile(delete=False) as f:
                test_file = Path(f.name)
                f.write(b"test content")

            try:
                test_data = (np.array([1, 2, 3]), np.array([4, 5, 6]), {})

                # Put should do nothing when disabled
                cache.put(test_file, test_data)
                assert len(cache.cache) == 0

                # Get should return None when disabled
                result = cache.get(test_file)
                assert result is None

            finally:
                test_file.unlink()

        finally:
            config.set("performance.cache_enabled", original_enabled)

    def test_cache_eviction(self):
        """Test cache eviction when size limit exceeded."""
        cache = DataCache(max_size_mb=0.001)  # Very small cache

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create multiple test files
            test_files = []
            test_data_sets = []

            for i in range(3):
                test_file = Path(temp_dir) / f"test_{i}.dat"
                test_file.write_text(f"test content {i}")
                test_files.append(test_file)

                # Large enough data to trigger eviction
                large_data = (
                    np.zeros(1000),  # Large array
                    np.zeros(1000),  # Large array
                    {"index": i},
                )
                test_data_sets.append(large_data)

            # Add data that should trigger evictions
            for file, data in zip(test_files, test_data_sets):
                cache.put(file, data)

            # Cache should be limited in size
            assert len(cache.cache) <= 2  # Some items should be evicted

    def test_cache_clear(self):
        """Test cache clearing."""
        cache = DataCache()

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test content")

        try:
            test_data = (np.array([1, 2, 3]), np.array([4, 5, 6]), {})

            # Add data to cache
            cache.put(test_file, test_data)
            assert len(cache.cache) == 1

            # Clear cache
            cache.clear()
            assert len(cache.cache) == 0
            assert len(cache.access_order) == 0

        finally:
            test_file.unlink()

    def test_cache_stats(self):
        """Test cache statistics."""
        cache = DataCache(max_size_mb=10)

        # Empty cache stats
        stats = cache.get_stats()
        assert stats["entries"] == 0
        assert stats["total_size_mb"] == 0
        assert stats["max_size_mb"] == 10
        assert stats["usage_percent"] == 0

        # Add some data
        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test content")

        try:
            test_data = (np.array([1, 2, 3]), np.array([4, 5, 6]), {})
            cache.put(test_file, test_data)

            stats = cache.get_stats()
            assert stats["entries"] == 1
            assert stats["total_size_mb"] > 0
            assert 0 < stats["usage_percent"] < 100

        finally:
            test_file.unlink()


class TestOptimizedLoader:
    """Test OptimizedLoader functionality."""

    def test_loader_initialization(self):
        """Test optimized loader initialization."""
        loader = OptimizedLoader(chunk_size_mb=5, cache_enabled=True)

        assert loader.chunk_size_mb == 5
        assert loader.cache is not None
        assert isinstance(loader.cache, DataCache)

    def test_loader_initialization_no_cache(self):
        """Test optimized loader without cache."""
        loader = OptimizedLoader(cache_enabled=False)

        assert loader.cache is None

    def test_loader_initialization_from_config(self):
        """Test loader initialization using config values."""
        from epyr.config import config

        original_chunk = config.get("performance.chunk_size_mb")
        try:
            config.set("performance.chunk_size_mb", 15)
            loader = OptimizedLoader()
            assert loader.chunk_size_mb == 15
        finally:
            config.set("performance.chunk_size_mb", original_chunk)

    @patch("epyr.performance.eprload")
    @patch("epyr.performance.MemoryMonitor.check_memory_limit")
    def test_load_epr_file_success(self, mock_memory_check, mock_eprload):
        """Test successful EPR file loading."""
        # Mock successful loading
        mock_memory_check.return_value = True
        mock_eprload.return_value = (
            np.linspace(3400, 3500, 100),
            np.random.randn(100),
            {"frequency": 9.4e9},
            "test.dsc",
        )

        loader = OptimizedLoader(cache_enabled=False)

        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test EPR data")

        try:
            result = loader.load_epr_file(test_file)

            assert result is not None
            assert len(result) == 4  # x, y, params, file_path
            mock_eprload.assert_called_once_with(str(test_file), plot_if_possible=False)

        finally:
            test_file.unlink()

    @patch("epyr.performance.eprload")
    def test_load_epr_file_with_cache(self, mock_eprload):
        """Test EPR file loading with cache."""
        test_data = (
            np.linspace(3400, 3500, 100),
            np.random.randn(100),
            {"frequency": 9.4e9},
            "test.dsc",
        )
        mock_eprload.return_value = test_data

        loader = OptimizedLoader(cache_enabled=True)

        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test EPR data")

        try:
            # First load - should call eprload and cache result
            result1 = loader.load_epr_file(test_file)
            assert mock_eprload.call_count == 1

            # Second load - should use cache
            result2 = loader.load_epr_file(test_file)
            assert mock_eprload.call_count == 1  # Not called again

            # Results should be identical
            np.testing.assert_array_equal(result1[0], result2[0])
            np.testing.assert_array_equal(result1[1], result2[1])

        finally:
            test_file.unlink()

    @patch("epyr.performance.eprload")
    @patch("epyr.performance.MemoryMonitor.check_memory_limit")
    @patch("epyr.performance.MemoryMonitor.optimize_memory")
    def test_load_epr_file_memory_optimization(
        self, mock_optimize, mock_memory_check, mock_eprload
    ):
        """Test memory optimization during loading."""
        # Mock memory limit exceeded
        mock_memory_check.return_value = False
        mock_eprload.return_value = (
            np.array([1, 2, 3]),
            np.array([4, 5, 6]),
            {},
            "test.dsc",
        )

        loader = OptimizedLoader()

        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test data")

        try:
            loader.load_epr_file(test_file)

            # Should have triggered memory optimization
            mock_optimize.assert_called_once()

        finally:
            test_file.unlink()

    @patch("epyr.performance.eprload")
    def test_load_epr_file_error_handling(self, mock_eprload):
        """Test error handling during file loading."""
        mock_eprload.side_effect = Exception("Loading failed")

        loader = OptimizedLoader()

        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = Path(f.name)
            f.write(b"test data")

        try:
            with pytest.raises(Exception):
                loader.load_epr_file(test_file)

        finally:
            test_file.unlink()

    def test_load_chunked_data_small_file(self):
        """Test chunked loading with small file."""
        loader = OptimizedLoader(chunk_size_mb=10)

        # Mock processor function
        mock_processor = Mock()
        mock_processor.return_value = "processed_result"

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"small file content")  # < 10 MB

        try:
            result = loader.load_chunked_data(test_file, mock_processor)

            # Should call processor with is_chunked=False for small file
            mock_processor.assert_called_once_with(test_file, is_chunked=False)
            assert result == "processed_result"

        finally:
            test_file.unlink()

    def test_load_chunked_data_large_file(self):
        """Test chunked loading with large file (warning case)."""
        loader = OptimizedLoader(chunk_size_mb=0.001)  # Very small chunk size

        # Mock processor function
        mock_processor = Mock()
        mock_processor.return_value = "processed_result"

        with tempfile.NamedTemporaryFile(delete=False) as f:
            test_file = Path(f.name)
            f.write(b"large file content that exceeds chunk size")

        try:
            # Should issue warning and fall back to non-chunked processing
            with pytest.warns(
                UserWarning, match="Chunked processing not yet implemented"
            ):
                result = loader.load_chunked_data(test_file, mock_processor)

            mock_processor.assert_called_once_with(test_file, is_chunked=False)
            assert result == "processed_result"

        finally:
            test_file.unlink()


class TestPerformanceUtilities:
    """Test performance utility functions."""

    @patch("epyr.performance.mkl")
    @patch("epyr.performance.os.cpu_count")
    def test_optimize_numpy_operations_with_mkl(self, mock_cpu_count, mock_mkl):
        """Test NumPy optimization with MKL available."""
        mock_cpu_count.return_value = 8

        from epyr.config import config

        original_parallel = config.get("performance.parallel_processing")

        try:
            config.set("performance.parallel_processing", True)

            optimize_numpy_operations()

            # Should set MKL threads to half of CPU count
            mock_mkl.set_num_threads.assert_called_with(4)

        finally:
            config.set("performance.parallel_processing", original_parallel)

    @patch("epyr.performance.os.cpu_count")
    def test_optimize_numpy_operations_without_mkl(self, mock_cpu_count):
        """Test NumPy optimization without MKL."""
        mock_cpu_count.return_value = 4

        from epyr.config import config

        original_parallel = config.get("performance.parallel_processing")

        try:
            config.set("performance.parallel_processing", True)

            # Mock MKL import failure
            with patch.dict("sys.modules", {"mkl": None}):
                optimize_numpy_operations()

            # Should set environment variable
            assert "OMP_NUM_THREADS" in os.environ

        finally:
            config.set("performance.parallel_processing", original_parallel)

    def test_optimize_numpy_operations_disabled(self):
        """Test NumPy optimization when parallel processing disabled."""
        from epyr.config import config

        original_parallel = config.get("performance.parallel_processing")

        try:
            config.set("performance.parallel_processing", False)

            # Should not modify anything when disabled
            optimize_numpy_operations()
            # No assertions needed - just ensure no exceptions

        finally:
            config.set("performance.parallel_processing", original_parallel)

    @patch("epyr.performance.MemoryMonitor.get_memory_info")
    @patch("epyr.performance.os.cpu_count")
    def test_get_performance_info(self, mock_cpu_count, mock_memory_info):
        """Test performance information gathering."""
        mock_cpu_count.return_value = 8
        mock_memory_info.return_value = {
            "rss_mb": 150.0,
            "vms_mb": 300.0,
            "percent": 3.5,
        }

        perf_info = get_performance_info()

        assert "memory" in perf_info
        assert "cpu_count" in perf_info
        assert "config" in perf_info

        assert perf_info["cpu_count"] == 8
        assert perf_info["memory"]["rss_mb"] == 150.0

        # Check config section
        config_info = perf_info["config"]
        assert "cache_enabled" in config_info
        assert "cache_size_mb" in config_info
        assert "chunk_size_mb" in config_info
        assert "memory_limit_mb" in config_info
        assert "parallel_processing" in config_info


class TestGlobalCacheInstance:
    """Test global cache instance."""

    def test_global_cache_import(self):
        """Test importing global cache instance."""
        from epyr.performance import _global_cache, get_global_cache

        assert _global_cache is not None
        assert isinstance(_global_cache, DataCache)

        # Should return same instance
        assert get_global_cache() is _global_cache
