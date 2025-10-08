#!/usr/bin/env python3
"""
EPyR Tools Demo 08: Performance Optimization and Caching
========================================================

This script demonstrates the performance optimization features of EPyR Tools
for handling large EPR datasets efficiently.

Functions demonstrated:
- MemoryMonitor - Memory usage monitoring and optimization
- DataCache - LRU caching for frequently accessed files
- OptimizedLoader - Memory-efficient data loading with chunking
- Performance benchmarking and profiling
"""

import sys
import time
import gc
from pathlib import Path
import numpy as np

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_memory_monitoring():
    """Demonstrate memory monitoring capabilities."""
    print("=== EPyR Tools Performance Demo - Memory Monitoring ===")
    print()

    print("1. Memory monitoring:")
    print("-" * 20)

    # Get initial memory info
    try:
        from epyr.performance import MemoryMonitor

        print("Initial memory status:")
        memory_info = MemoryMonitor.get_memory_info()
        print(f"  RSS (Resident Set Size): {memory_info['rss_mb']:.1f} MB")
        print(f"  VMS (Virtual Memory Size): {memory_info['vms_mb']:.1f} MB")
        print(f"  Memory percentage: {memory_info['percent']:.1f}%")

        # Check memory limit
        within_limit = MemoryMonitor.check_memory_limit()
        print(f"  Within configured memory limit: {within_limit}")

        print()

        # Simulate memory usage with large arrays
        print("Simulating memory usage with large arrays:")
        large_arrays = []

        for i in range(3):
            # Create 10MB array
            array_size = 10 * 1024 * 1024 // 8  # 10 MB in float64 elements
            large_array = np.random.random(array_size)
            large_arrays.append(large_array)

            memory_info = MemoryMonitor.get_memory_info()
            print(f"  After array {i+1}: {memory_info['rss_mb']:.1f} MB RSS")

        # Optimize memory
        print("\nOptimizing memory:")
        MemoryMonitor.optimize_memory()
        del large_arrays
        gc.collect()

        memory_info = MemoryMonitor.get_memory_info()
        print(f"  After optimization: {memory_info['rss_mb']:.1f} MB RSS")

    except ImportError as e:
        print(f"Performance monitoring not available: {e}")
        print("Install psutil for full memory monitoring capabilities")

    print()


def demo_data_caching():
    """Demonstrate data caching system."""
    print("2. Data caching system:")
    print("-" * 23)

    try:
        from epyr.performance import DataCache

        # Initialize cache with small size for demo
        cache = DataCache(max_size_mb=50)
        print(f"Initialized cache with {cache.max_size_mb} MB limit")

        # Create synthetic data files for caching demo
        cache_dir = Path(__file__).parent / "cache_demo"
        cache_dir.mkdir(exist_ok=True)

        synthetic_files = []
        for i in range(3):
            filename = f"synthetic_data_{i+1}.npy"
            filepath = cache_dir / filename

            # Create different sized arrays
            size = (i + 1) * 1000
            data = np.random.random(size)
            np.save(filepath, data)
            synthetic_files.append(filepath)
            print(f"  Created {filename}: {data.nbytes/1024:.1f} KB")

        print()

        # Demonstrate caching with load function
        def load_numpy_file(filepath):
            """Simple load function for demonstration."""
            print(f"    Loading {filepath.name} from disk...")
            return np.load(filepath)

        # Load files and cache them
        print("Loading files with caching:")
        cache_times = []
        disk_times = []

        for filepath in synthetic_files:
            # First load (from disk)
            start_time = time.time()

            # Check cache first
            cached_data = cache.get(filepath)
            if cached_data is None:
                # Load from disk and cache
                data = load_numpy_file(filepath)
                cache.put(filepath, (data,))
            else:
                data = cached_data[0]

            disk_time = time.time() - start_time
            disk_times.append(disk_time)
            print(f"  First load {filepath.name}: {disk_time*1000:.2f} ms")

            # Second load (from cache)
            start_time = time.time()
            cached_data = cache.get(filepath)
            if cached_data is not None:
                cached_data = cached_data[0]
            else:
                cached_data = load_numpy_file(filepath)
                cache.put(filepath, (cached_data,))

            cache_time = time.time() - start_time
            cache_times.append(cache_time)
            print(f"  Cached load {filepath.name}: {cache_time*1000:.2f} ms")

            # Verify data integrity
            assert np.array_equal(data, cached_data), "Cache data mismatch!"

        print()
        print("Cache performance summary:")
        avg_disk_time = np.mean(disk_times) * 1000
        avg_cache_time = np.mean(cache_times) * 1000
        speedup = avg_disk_time / avg_cache_time if avg_cache_time > 0 else float('inf')
        print(f"  Average disk load time: {avg_disk_time:.2f} ms")
        print(f"  Average cache load time: {avg_cache_time:.2f} ms")
        print(f"  Cache speedup: {speedup:.1f}x")

        # Show cache statistics
        print(f"\nCache statistics:")
        cache_size = sum(entry['size_mb'] for entry in cache.cache.values()) if hasattr(cache, 'cache') else 0
        print(f"  Cache entries: {len(cache.cache) if hasattr(cache, 'cache') else 0}")
        print(f"  Max cache size: {cache.max_size_mb} MB")
        print(f"  Current cache size: {cache_size:.1f} MB")

        # Clean up
        for filepath in synthetic_files:
            filepath.unlink()
        cache_dir.rmdir()

    except ImportError as e:
        print(f"Caching system not available: {e}")

    print()


def demo_optimized_loader():
    """Demonstrate optimized data loading."""
    print("3. Optimized data loading:")
    print("-" * 26)

    try:
        from epyr.performance import OptimizedLoader

        # Initialize optimized loader
        loader = OptimizedLoader(chunk_size_mb=5, cache_enabled=True)
        print("Initialized OptimizedLoader with 5MB chunks and caching enabled")

        data_dir = Path(__file__).parent.parent.parent / "data"

        # Look for real EPR files
        epr_files = []
        for pattern in ["*.DSC", "*.dsc", "*.PAR", "*.par"]:
            epr_files.extend(data_dir.glob(pattern))

        if epr_files:
            test_file = epr_files[0]
            print(f"Testing with real file: {test_file.name}")

            # Time optimized loading
            start_time = time.time()
            try:
                x, y, params = loader.load_epr_file(str(test_file))
                load_time = time.time() - start_time

                print(f"  Load time: {load_time*1000:.1f} ms")
                print(f"  Data shape: {y.shape}")
                print(f"  Memory usage optimized: {loader.memory_optimized}")

                # Load again to test caching
                start_time = time.time()
                x2, y2, params2 = loader.load_epr_file(str(test_file))
                cached_load_time = time.time() - start_time

                print(f"  Cached load time: {cached_load_time*1000:.1f} ms")
                print(f"  Cache speedup: {load_time/cached_load_time:.1f}x")

            except Exception as e:
                print(f"  Error with optimized loader: {e}")

        else:
            print("No real EPR files found, creating synthetic large dataset...")
            create_synthetic_performance_test(loader)

    except ImportError as e:
        print(f"OptimizedLoader not available: {e}")

    print()


def create_synthetic_performance_test(loader=None):
    """Create synthetic large dataset for performance testing."""
    print("Creating synthetic large dataset for performance testing:")

    # Create large 2D dataset
    print("  Generating large 2D EPR dataset (simulated angular study)...")
    n_angles = 180
    n_field_points = 1024
    total_points = n_angles * n_field_points

    print(f"  Dataset size: {n_angles} x {n_field_points} = {total_points:,} points")

    # Simulate memory usage
    estimated_memory_mb = total_points * 8 / (1024 * 1024)  # float64
    print(f"  Estimated memory: {estimated_memory_mb:.1f} MB")

    # Generate data in chunks to simulate large file loading
    chunk_size = 10000  # Points per chunk
    n_chunks = total_points // chunk_size + 1

    print(f"  Processing in {n_chunks} chunks of {chunk_size} points each")

    start_time = time.time()
    processed_points = 0

    for i in range(n_chunks):
        chunk_start = i * chunk_size
        chunk_end = min(chunk_start + chunk_size, total_points)
        chunk_points = chunk_end - chunk_start

        if chunk_points <= 0:
            break

        # Simulate processing
        chunk_data = np.random.random(chunk_points)
        # Simulate some computation
        processed_chunk = np.sin(chunk_data) * np.exp(-chunk_data)
        processed_points += chunk_points

        if i % 10 == 0:  # Progress update every 10 chunks
            progress = processed_points / total_points * 100
            elapsed = time.time() - start_time
            print(f"    Progress: {progress:.1f}% ({elapsed:.1f}s)")

    total_time = time.time() - start_time
    points_per_second = total_points / total_time
    print(f"  Total processing time: {total_time:.2f} seconds")
    print(f"  Processing rate: {points_per_second:,.0f} points/second")


def demo_performance_benchmarking():
    """Demonstrate performance benchmarking."""
    print("4. Performance benchmarking:")
    print("-" * 28)

    # Benchmark different data operations
    operations = {
        'Array creation': lambda n: np.random.random(n),
        'FFT': lambda n: np.fft.fft(np.random.random(n)),
        'Matrix multiplication': lambda n: np.dot(np.random.random((int(np.sqrt(n)), int(np.sqrt(n)))),
                                                 np.random.random((int(np.sqrt(n)), int(np.sqrt(n))))),
        'Polynomial fit': lambda n: np.polyfit(np.arange(n), np.random.random(n), 3)
    }

    sizes = [1000, 10000, 100000]

    print("Benchmarking common EPR data operations:")
    print(f"{'Operation':<20} {'1K pts':<10} {'10K pts':<10} {'100K pts':<10}")
    print("-" * 55)

    for op_name, op_func in operations.items():
        times = []
        for size in sizes:
            try:
                # Adjust size for matrix operations
                if 'Matrix' in op_name and size > 10000:
                    continue  # Skip large matrix operations

                start_time = time.time()
                result = op_func(size)
                op_time = time.time() - start_time
                times.append(f"{op_time*1000:.1f}ms")
            except Exception as e:
                times.append("N/A")

        # Fill missing times
        while len(times) < 3:
            times.append("N/A")

        print(f"{op_name:<20} {times[0]:<10} {times[1]:<10} {times[2]:<10}")

    print()


def demo_configuration_performance():
    """Demonstrate performance configuration options."""
    print("5. Performance configuration:")
    print("-" * 29)

    print("Current performance settings:")

    # Show performance-related configuration
    performance_keys = [
        'performance.cache_enabled',
        'performance.cache_size_mb',
        'performance.memory_limit_mb',
        'performance.chunk_size_mb',
        'performance.parallel_processing'
    ]

    for key in performance_keys:
        try:
            value = epyr.config.get(key, "Not configured")
            print(f"  {key}: {value}")
        except:
            print(f"  {key}: Configuration not accessible")

    print()

    # Demonstrate configuration changes
    print("Demonstrating configuration changes:")
    try:
        # Show current cache setting
        current_cache = epyr.config.get('performance.cache_enabled', True)
        print(f"  Cache currently: {'enabled' if current_cache else 'disabled'}")

        # Temporarily change setting
        epyr.config.set('performance.cache_enabled', not current_cache)
        new_cache = epyr.config.get('performance.cache_enabled')
        print(f"  Cache changed to: {'enabled' if new_cache else 'disabled'}")

        # Restore original setting
        epyr.config.set('performance.cache_enabled', current_cache)
        print(f"  Cache restored to: {'enabled' if current_cache else 'disabled'}")

    except Exception as e:
        print(f"  Configuration change demo failed: {e}")

    print()


def demo_memory_profiling():
    """Demonstrate memory profiling of EPR operations."""
    print("6. Memory profiling of EPR operations:")
    print("-" * 38)

    try:
        from epyr.performance import MemoryMonitor

        # Profile different EPR operations
        operations = [
            ("Load synthetic data", lambda: np.random.random(100000)),
            ("Calculate FFT", lambda: np.fft.fft(np.random.random(50000))),
            ("Baseline fit", lambda: np.polyfit(np.arange(10000), np.random.random(10000), 5)),
            ("Gaussian convolution", lambda: np.convolve(np.random.random(10000),
                                                       np.exp(-np.arange(-50, 51)**2/100), mode='same'))
        ]

        print("Memory usage during common EPR operations:")
        print(f"{'Operation':<20} {'Before (MB)':<12} {'After (MB)':<12} {'Change (MB)':<12}")
        print("-" * 60)

        for op_name, op_func in operations:
            # Get memory before operation
            mem_before = MemoryMonitor.get_memory_info()

            # Perform operation
            result = op_func()

            # Get memory after operation
            mem_after = MemoryMonitor.get_memory_info()

            change = mem_after['rss_mb'] - mem_before['rss_mb']

            print(f"{op_name:<20} {mem_before['rss_mb']:<12.1f} {mem_after['rss_mb']:<12.1f} {change:<12.1f}")

            # Clean up
            del result
            gc.collect()

    except ImportError:
        print("Memory profiling requires psutil package")

    print()


def main():
    """Run all performance demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_memory_monitoring()
    demo_data_caching()
    demo_optimized_loader()
    demo_performance_benchmarking()
    demo_configuration_performance()
    demo_memory_profiling()

    print("=== Performance Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- MemoryMonitor provides real-time memory usage tracking and optimization")
    print("- DataCache implements LRU caching for frequently accessed files")
    print("- OptimizedLoader handles large datasets with chunking and memory management")
    print("- Performance can be tuned through configuration settings")
    print("- Memory profiling helps identify resource-intensive operations")
    print("- Caching can provide significant speedup for repeated file access")
    print("- Proper memory management is crucial for large EPR datasets")
    print()
    print("Performance optimization tips:")
    print("- Enable caching for workflows with repeated file access")
    print("- Use chunked processing for files larger than available memory")
    print("- Monitor memory usage to prevent system overload")
    print("- Configure cache size based on available system memory")


if __name__ == "__main__":
    main()