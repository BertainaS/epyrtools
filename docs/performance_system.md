# Performance System - epyr.performance

The `epyr.performance` module provides comprehensive performance optimization tools for handling large EPR datasets efficiently, including memory management, caching, and optimized data loading.

## Overview

EPyR Tools includes sophisticated performance optimization features:
- **Memory Management**: Intelligent monitoring and optimization
- **Data Caching**: LRU cache for frequently accessed files
- **Optimized Loading**: Chunked processing for large datasets
- **System Monitoring**: Resource usage tracking and alerts
- **NumPy Optimization**: Multi-core processing configuration

## Architecture

### Core Components

1. **MemoryMonitor**: System memory tracking and optimization
2. **DataCache**: Intelligent caching for EPR data files
3. **OptimizedLoader**: Enhanced data loading with performance features
4. **Performance Utilities**: System optimization functions

## Memory Management

### MemoryMonitor Class

Monitors and optimizes memory usage during data processing:

```python
from epyr.performance import MemoryMonitor

# Get current memory usage
memory_info = MemoryMonitor.get_memory_info()
print(f"Memory usage: {memory_info['rss_mb']:.1f} MB")
print(f"Virtual memory: {memory_info['vms_mb']:.1f} MB")
print(f"Memory percentage: {memory_info['percent']:.1f}%")

# Check if within limits
if not MemoryMonitor.check_memory_limit():
    print("Memory limit exceeded!")
    MemoryMonitor.optimize_memory()  # Force garbage collection
```

### Memory Information Structure

```python
{
    'rss_mb': 150.5,      # Resident Set Size in MB
    'vms_mb': 300.2,      # Virtual Memory Size in MB  
    'percent': 3.8        # Percentage of system memory used
}
```

### Memory Limit Checking

```python
# Automatic memory monitoring
def process_large_dataset():
    if not MemoryMonitor.check_memory_limit():
        print("Approaching memory limit, optimizing...")
        MemoryMonitor.optimize_memory()
        
        # Still over limit? Reduce processing parameters
        if not MemoryMonitor.check_memory_limit():
            chunk_size = chunk_size // 2  # Reduce chunk size
```

### Memory Optimization

The system automatically:
- Runs garbage collection when limits are approached
- Provides warnings when memory usage is high
- Integrates with configuration system for limit settings

```python
from epyr.config import config

# Configure memory limits
config.set('performance.memory_limit_mb', 1000)  # 1GB limit
```

## Data Caching System

### DataCache Class

Intelligent LRU (Least Recently Used) cache for EPR data files:

```python
from epyr.performance import DataCache

# Create cache with 200MB limit
cache = DataCache(max_size_mb=200)

# Cache data automatically considers file modification times
data = (x_array, y_array, parameters_dict)
cache.put(file_path, data)

# Retrieve cached data (None if not found or invalid)
cached_data = cache.get(file_path)
if cached_data:
    x, y, params = cached_data
    print("Using cached data")
else:
    print("Loading fresh data")
```

### Cache Features

**Automatic File Tracking:**
```python
# Cache keys include file modification time and size
# Data is automatically invalidated when file changes
cache_key = file_path + "_" + modification_time + "_" + file_size
```

**Memory Management:**
```python
# Automatic eviction when cache limit reached
cache = DataCache(max_size_mb=100)

# Add large dataset
large_data = (np.zeros(1000000), np.zeros(1000000), {})
cache.put(Path("large_file.dsc"), large_data)

# Older entries automatically evicted to make space
stats = cache.get_stats()
print(f"Cache usage: {stats['usage_percent']:.1f}%")
```

**Cache Statistics:**
```python
stats = cache.get_stats()
# Returns:
{
    'entries': 5,                # Number of cached files
    'total_size_mb': 45.2,      # Total cache size
    'max_size_mb': 100,         # Cache limit
    'usage_percent': 45.2       # Percentage used
}
```

### Global Cache Instance

```python
from epyr.performance import get_global_cache

# Access the global cache instance
cache = get_global_cache()

# All EPyR Tools modules can use this shared cache
cache.clear()  # Clear all cached data
```

## Optimized Data Loading

### OptimizedLoader Class

Enhanced EPR data loader with performance optimizations:

```python
from epyr.performance import OptimizedLoader

# Create optimized loader
loader = OptimizedLoader(
    chunk_size_mb=10,      # Process in 10MB chunks
    cache_enabled=True     # Enable caching
)

# Load EPR file with optimizations
x, y, params, file_path = loader.load_epr_file('spectrum.dsc')
```

### Key Features

**Memory-Efficient Loading:**
```python
# Automatic memory checking before loading
loader = OptimizedLoader()

# Memory optimization is automatic
data = loader.load_epr_file('large_spectrum.dsc')
# - Checks memory before loading
# - Optimizes memory if limits exceeded  
# - Uses chunked processing for very large files
```

**Intelligent Caching:**
```python
# Second load of same file uses cache
data1 = loader.load_epr_file('spectrum.dsc')  # Loads from disk
data2 = loader.load_epr_file('spectrum.dsc')  # Uses cache

# Cache automatically detects file changes
# Modifying file invalidates cache entry
```

**Chunked Processing Framework:**
```python
def custom_processor(file_path, is_chunked=False):
    """Custom processing function for chunked data."""
    if is_chunked:
        return process_in_chunks(file_path)
    else:
        return process_entire_file(file_path)

# Use chunked processing for large files
result = loader.load_chunked_data('huge_file.dsc', custom_processor)
```

### Configuration Integration

```python
from epyr.config import config

# Loader uses configuration automatically
config.set('performance.chunk_size_mb', 20)
config.set('performance.cache_enabled', True)

loader = OptimizedLoader()  # Uses configured settings
```

## System Optimization

### NumPy Performance

```python
from epyr.performance import optimize_numpy_operations

# Optimize NumPy for multi-core processing
optimize_numpy_operations()
```

**Automatic Optimizations:**
- Sets optimal thread count for MKL operations
- Configures OpenMP threads for NumPy
- Uses configuration settings for parallel processing control

```python
# Configuration controls optimization
config.set('performance.parallel_processing', True)
optimize_numpy_operations()  # Uses half of available cores
```

### Performance Information

```python
from epyr.performance import get_performance_info

perf_info = get_performance_info()
```

**Information Includes:**
```python
{
    'memory': {
        'rss_mb': 150.5,
        'vms_mb': 300.2, 
        'percent': 3.8
    },
    'cpu_count': 8,
    'config': {
        'cache_enabled': True,
        'cache_size_mb': 100,
        'chunk_size_mb': 10,
        'memory_limit_mb': 500,
        'parallel_processing': True
    }
}
```

## Integration Examples

### Module Integration

```python
# In any EPyR Tools module
from epyr.performance import OptimizedLoader, MemoryMonitor
from epyr.config import config

def load_epr_data_optimized(file_path):
    """Load EPR data with full optimization."""
    
    # Check memory before starting
    if not MemoryMonitor.check_memory_limit():
        logger.warning("Memory usage high, optimizing...")
        MemoryMonitor.optimize_memory()
    
    # Use optimized loader
    loader = OptimizedLoader(
        chunk_size_mb=config.get('performance.chunk_size_mb'),
        cache_enabled=config.get('performance.cache_enabled')
    )
    
    return loader.load_epr_file(file_path)
```

### CLI Integration

```python
# In CLI commands
def cmd_batch_convert():
    # Use performance monitoring for batch operations
    total_files = len(files_to_process)
    
    for i, file_path in enumerate(files_to_process):
        # Monitor memory during batch processing
        if i % 10 == 0:  # Check every 10 files
            if not MemoryMonitor.check_memory_limit():
                MemoryMonitor.optimize_memory()
        
        # Process with optimization
        loader = OptimizedLoader()
        data = loader.load_epr_file(file_path)
        process_data(data)
```

### FAIR Conversion Integration

```python
# In FAIR conversion module
from epyr.performance import OptimizedLoader

def convert_bruker_to_fair(input_file, **kwargs):
    """FAIR conversion with performance optimization."""
    
    # Use optimized loading for better performance
    loader = OptimizedLoader(cache_enabled=True)
    
    try:
        x, y, params, file_path = loader.load_epr_file(input_file)
        # Conversion process...
        
    except MemoryError:
        # Handle memory issues gracefully
        logger.warning("Memory limit reached, trying with smaller chunks")
        loader = OptimizedLoader(chunk_size_mb=5)  # Smaller chunks
        x, y, params, file_path = loader.load_epr_file(input_file)
```

## Advanced Usage

### Custom Cache Management

```python
from epyr.performance import DataCache

class CustomDataCache(DataCache):
    """Custom cache with additional features."""
    
    def __init__(self, max_size_mb=100):
        super().__init__(max_size_mb)
        self.hit_count = 0
        self.miss_count = 0
    
    def get(self, file_path):
        result = super().get(file_path)
        if result:
            self.hit_count += 1
        else:
            self.miss_count += 1
        return result
    
    def get_hit_rate(self):
        total = self.hit_count + self.miss_count
        return self.hit_count / total if total > 0 else 0

# Use custom cache
custom_cache = CustomDataCache(max_size_mb=200)
```

### Performance Monitoring

```python
import time
from epyr.performance import MemoryMonitor, get_performance_info

def monitor_operation(operation_name, func, *args, **kwargs):
    """Monitor performance of an operation."""
    
    # Get initial state
    start_time = time.time()
    start_memory = MemoryMonitor.get_memory_info()
    
    print(f"Starting {operation_name}")
    print(f"Initial memory: {start_memory['rss_mb']:.1f} MB")
    
    # Run operation
    result = func(*args, **kwargs)
    
    # Get final state
    end_time = time.time()
    end_memory = MemoryMonitor.get_memory_info()
    
    # Report results
    duration = end_time - start_time
    memory_delta = end_memory['rss_mb'] - start_memory['rss_mb']
    
    print(f"Operation completed in {duration:.2f}s")
    print(f"Memory change: {memory_delta:+.1f} MB")
    print(f"Final memory: {end_memory['rss_mb']:.1f} MB")
    
    return result

# Use monitoring
result = monitor_operation("EPR Loading", load_epr_data, "spectrum.dsc")
```

### Large Dataset Processing

```python
from epyr.performance import OptimizedLoader, MemoryMonitor

def process_large_dataset(file_paths):
    """Process multiple large EPR files efficiently."""
    
    loader = OptimizedLoader(
        chunk_size_mb=5,  # Small chunks for large files
        cache_enabled=True
    )
    
    results = []
    for i, file_path in enumerate(file_paths):
        print(f"Processing {i+1}/{len(file_paths)}: {file_path.name}")
        
        # Check memory every few files
        if i % 5 == 0:
            memory_info = MemoryMonitor.get_memory_info()
            print(f"Memory usage: {memory_info['rss_mb']:.1f} MB")
            
            if not MemoryMonitor.check_memory_limit():
                print("Optimizing memory...")
                MemoryMonitor.optimize_memory()
                
                # Clear cache if still over limit
                if not MemoryMonitor.check_memory_limit():
                    loader.cache.clear()
                    print("Cache cleared due to memory pressure")
        
        # Process file
        try:
            data = loader.load_epr_file(file_path)
            result = analyze_epr_data(data)
            results.append(result)
            
        except MemoryError:
            print(f"Skipping {file_path} due to memory constraints")
            continue
    
    return results
```

## Configuration Options

### Performance Settings

All performance features are configurable:

```python
from epyr.config import config

# Cache configuration
config.set('performance.cache_enabled', True)
config.set('performance.cache_size_mb', 200)

# Memory management
config.set('performance.memory_limit_mb', 1000)
config.set('performance.chunk_size_mb', 20)

# Parallel processing
config.set('performance.parallel_processing', True)

# Save configuration
config.save()
```

### Environment Variables

Performance can be controlled via environment variables:

```bash
export EPYR_PERFORMANCE_CACHE_ENABLED=true
export EPYR_PERFORMANCE_CACHE_SIZE_MB=200
export EPYR_PERFORMANCE_MEMORY_LIMIT_MB=1000
export EPYR_PERFORMANCE_CHUNK_SIZE_MB=20
export EPYR_PERFORMANCE_PARALLEL_PROCESSING=true
```

## Best Practices

### Memory Management
- Monitor memory usage for large batch operations
- Use appropriate cache sizes based on available system memory
- Clear cache when processing many different files
- Implement graceful degradation when memory limits are reached

### Caching Strategy
- Enable caching for workflows that reload the same files
- Use larger cache sizes for analysis workflows
- Disable caching for one-time conversion operations
- Monitor cache hit rates to optimize cache size

### Large Dataset Handling
- Use chunked processing for files > 100MB
- Reduce chunk sizes when memory is constrained
- Process files in order of increasing size when possible
- Implement progress reporting for long operations

### Performance Monitoring
- Use performance info to diagnose bottlenecks
- Monitor memory usage during development
- Profile operations to optimize critical paths
- Log performance metrics for production deployments

The performance system provides comprehensive tools for handling EPR data efficiently, from small single files to large batch processing operations, while maintaining system stability and optimal resource usage.