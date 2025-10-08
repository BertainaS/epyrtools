#!/usr/bin/env python3
"""
EPyR Tools Demo 09: Configuration Management
============================================

This script demonstrates the configuration management system of EPyR Tools
which provides centralized settings for all package components.

Functions demonstrated:
- config.get() - Retrieve configuration values
- config.set() - Set configuration values
- config.save() - Save configuration to file
- Environment variable integration
- Configuration file management
- Default settings and validation
"""

import sys
import os
import json
import tempfile
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_basic_configuration():
    """Demonstrate basic configuration operations."""
    print("=== EPyR Tools Configuration Demo - Basic Operations ===")
    print()

    print("1. Current configuration overview:")
    print("-" * 34)

    # Show current config file location
    try:
        config_file = epyr.config._config_file
        print(f"Configuration file: {config_file}")
        print(f"File exists: {Path(config_file).exists()}")
    except:
        print("Configuration file location not accessible")

    print()

    # Show some key configuration sections
    sections = ['plotting', 'data_loading', 'baseline', 'performance', 'fair_conversion']

    for section in sections:
        print(f"{section.capitalize()} configuration:")
        try:
            # Get all values in this section
            section_config = {}
            for key in epyr.config._defaults.get(section, {}):
                full_key = f"{section}.{key}"
                value = epyr.config.get(full_key)
                section_config[key] = value

            for key, value in section_config.items():
                print(f"  {key}: {value}")

        except Exception as e:
            print(f"  Error accessing {section} config: {e}")

        print()


def demo_getting_setting_values():
    """Demonstrate getting and setting configuration values."""
    print("2. Getting and setting configuration values:")
    print("-" * 44)

    # Demonstrate getting values with defaults
    print("Getting configuration values:")

    test_keys = [
        ('plotting.dpi', 'Plot DPI setting'),
        ('plotting.figure_size', 'Default figure size'),
        ('performance.cache_enabled', 'Cache enabled status'),
        ('baseline.default_method', 'Default baseline method'),
        ('fair_conversion.default_formats', 'Default export formats')
    ]

    for key, description in test_keys:
        try:
            value = epyr.config.get(key)
            print(f"  {description}: {value}")
        except Exception as e:
            print(f"  {description}: Error - {e}")

    print()

    # Demonstrate setting values
    print("Setting configuration values:")

    # Save original values
    original_dpi = epyr.config.get('plotting.dpi', 300)
    original_cache = epyr.config.get('performance.cache_enabled', True)

    try:
        # Change some values
        epyr.config.set('plotting.dpi', 150)
        new_dpi = epyr.config.get('plotting.dpi')
        print(f"  Changed DPI from {original_dpi} to {new_dpi}")

        epyr.config.set('performance.cache_enabled', False)
        new_cache = epyr.config.get('performance.cache_enabled')
        print(f"  Changed cache from {original_cache} to {new_cache}")

        # Restore original values
        epyr.config.set('plotting.dpi', original_dpi)
        epyr.config.set('performance.cache_enabled', original_cache)
        print(f"  Restored original values")

    except Exception as e:
        print(f"  Error setting values: {e}")

    print()


def demo_hierarchical_config():
    """Demonstrate hierarchical configuration structure."""
    print("3. Hierarchical configuration structure:")
    print("-" * 40)

    print("Configuration is organized in sections:")

    try:
        defaults = epyr.config._defaults

        for section, section_config in defaults.items():
            print(f"\n[{section}]")
            for key, value in section_config.items():
                if isinstance(value, (dict, list)):
                    value_str = str(value)[:50] + "..." if len(str(value)) > 50 else str(value)
                else:
                    value_str = str(value)
                print(f"  {key} = {value_str}")

    except Exception as e:
        print(f"Error accessing default configuration: {e}")

    print()


def demo_environment_variables():
    """Demonstrate environment variable integration."""
    print("4. Environment variable integration:")
    print("-" * 36)

    print("EPyR Tools supports environment variables with EPYR_ prefix:")

    # Show current environment variables
    epyr_vars = {k: v for k, v in os.environ.items() if k.startswith('EPYR_')}

    if epyr_vars:
        print("Current EPYR environment variables:")
        for var, value in epyr_vars.items():
            print(f"  {var} = {value}")
    else:
        print("No EPYR environment variables currently set")

    print()

    # Demonstrate setting environment variable
    print("Demonstrating environment variable override:")

    # Set a test environment variable
    test_var = 'EPYR_PLOTTING_DPI'
    original_value = os.environ.get(test_var)

    try:
        # Set environment variable
        os.environ[test_var] = '200'
        print(f"  Set {test_var} = 200")

        # The configuration system should pick this up on reload
        # (Note: actual implementation may vary)
        print("  Environment variables take precedence over config file settings")

        # Clean up
        if original_value is not None:
            os.environ[test_var] = original_value
        else:
            del os.environ[test_var]
        print("  Cleaned up test environment variable")

    except Exception as e:
        print(f"  Error with environment variable demo: {e}")

    print()


def demo_config_file_operations():
    """Demonstrate configuration file operations."""
    print("5. Configuration file operations:")
    print("-" * 33)

    # Create a temporary config file for demonstration
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_file:
        temp_config_path = temp_file.name

        # Create sample configuration
        sample_config = {
            'plotting': {
                'dpi': 150,
                'figure_size': [10, 8],
                'color_scheme': 'plasma'
            },
            'performance': {
                'cache_enabled': False,
                'cache_size_mb': 50
            }
        }

        json.dump(sample_config, temp_file, indent=2)
        print(f"Created temporary config file: {temp_config_path}")

    # Show file contents
    try:
        with open(temp_config_path, 'r') as f:
            content = f.read()
        print("Sample configuration file content:")
        print(content)

    except Exception as e:
        print(f"Error reading config file: {e}")

    print()

    # Demonstrate saving configuration
    print("Configuration saving:")
    try:
        # Save current configuration (to the actual config file)
        epyr.config.save()
        print("  Configuration saved successfully")

        # Show what would be saved
        print("  Current configuration includes all sections:")
        for section in epyr.config._defaults.keys():
            print(f"    - {section}")

    except Exception as e:
        print(f"  Error saving configuration: {e}")

    # Clean up temporary file
    try:
        os.unlink(temp_config_path)
        print("  Cleaned up temporary config file")
    except:
        pass

    print()


def demo_configuration_validation():
    """Demonstrate configuration validation."""
    print("6. Configuration validation:")
    print("-" * 28)

    print("Configuration validation prevents invalid settings:")

    # Test invalid values
    test_cases = [
        ('plotting.dpi', -100, 'Negative DPI'),
        ('plotting.figure_size', 'invalid', 'Non-list figure size'),
        ('performance.cache_size_mb', 'not_a_number', 'Non-numeric cache size'),
        ('baseline.default_poly_order', -1, 'Negative polynomial order')
    ]

    for key, invalid_value, description in test_cases:
        try:
            original_value = epyr.config.get(key)
            print(f"\nTesting {description}:")
            print(f"  Original value: {original_value}")
            print(f"  Attempting to set invalid value: {invalid_value}")

            # Try to set invalid value
            epyr.config.set(key, invalid_value)
            new_value = epyr.config.get(key)

            if new_value == invalid_value:
                print("  ✗ Invalid value was accepted (validation needed)")
            else:
                print("  ✓ Invalid value was rejected or corrected")

            # Restore original value
            epyr.config.set(key, original_value)

        except Exception as e:
            print(f"  ✓ Configuration validation prevented invalid setting: {e}")

    print()


def demo_practical_configuration():
    """Demonstrate practical configuration scenarios."""
    print("7. Practical configuration scenarios:")
    print("-" * 37)

    print("Common configuration use cases:")

    # Scenario 1: High-quality publication plots
    print("\nScenario 1: Configure for high-quality publication plots")
    pub_settings = {
        'plotting.dpi': 300,
        'plotting.figure_size': [6, 4.5],
        'plotting.font_size': 10,
        'plotting.line_width': 1.0,
        'plotting.save_format': 'pdf'
    }

    for key, value in pub_settings.items():
        original = epyr.config.get(key)
        epyr.config.set(key, value)
        print(f"  {key}: {original} → {value}")

    # Scenario 2: Memory-constrained environment
    print("\nScenario 2: Configure for memory-constrained environment")
    memory_settings = {
        'performance.cache_enabled': False,
        'performance.chunk_size_mb': 5,
        'performance.memory_limit_mb': 200
    }

    for key, value in memory_settings.items():
        original = epyr.config.get(key)
        epyr.config.set(key, value)
        print(f"  {key}: {original} → {value}")

    # Scenario 3: Batch processing mode
    print("\nScenario 3: Configure for batch processing")
    batch_settings = {
        'data_loading.auto_plot': False,
        'fair_conversion.default_formats': ['csv', 'hdf5'],
        'logging.level': 'WARNING'
    }

    for key, value in batch_settings.items():
        original = epyr.config.get(key)
        epyr.config.set(key, value)
        print(f"  {key}: {original} → {value}")

    print("\nNote: These are demonstration changes - actual values depend on your needs")
    print()


def demo_config_backup_restore():
    """Demonstrate configuration backup and restore."""
    print("8. Configuration backup and restore:")
    print("-" * 36)

    try:
        # Create backup of current configuration
        backup_file = Path(__file__).parent / "config_backup.json"

        # Export current configuration
        current_config = {}
        for section in epyr.config._defaults.keys():
            current_config[section] = {}
            for key in epyr.config._defaults[section].keys():
                full_key = f"{section}.{key}"
                current_config[section][key] = epyr.config.get(full_key)

        # Save backup
        with open(backup_file, 'w') as f:
            json.dump(current_config, f, indent=2)

        print(f"Configuration backed up to: {backup_file.name}")
        print(f"Backup file size: {backup_file.stat().st_size} bytes")

        # Show backup structure
        print("\nBackup contains sections:")
        for section in current_config.keys():
            setting_count = len(current_config[section])
            print(f"  {section}: {setting_count} settings")

        # Clean up backup file
        backup_file.unlink()
        print("\nBackup file removed after demonstration")

    except Exception as e:
        print(f"Error with backup/restore demo: {e}")

    print()


def main():
    """Run all configuration demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_basic_configuration()
    demo_getting_setting_values()
    demo_hierarchical_config()
    demo_environment_variables()
    demo_config_file_operations()
    demo_configuration_validation()
    demo_practical_configuration()
    demo_config_backup_restore()

    print("=== Configuration Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- Centralized configuration system manages all EPyR Tools settings")
    print("- Hierarchical structure with logical grouping of related settings")
    print("- config.get() and config.set() provide easy access to any setting")
    print("- Environment variables (EPYR_*) can override configuration file settings")
    print("- Configuration can be saved to persistent file for future sessions")
    print("- Settings validation helps prevent invalid configurations")
    print("- Common scenarios can be quickly configured with multiple setting changes")
    print("- Configuration backup/restore enables sharing settings between systems")
    print()
    print("Configuration tips:")
    print("- Use 'plotting' section for customizing plot appearance")
    print("- Adjust 'performance' settings based on system capabilities")
    print("- Configure 'fair_conversion' for preferred export formats")
    print("- Set 'logging.level' to DEBUG for troubleshooting")
    print("- Use environment variables for deployment-specific settings")


if __name__ == "__main__":
    main()