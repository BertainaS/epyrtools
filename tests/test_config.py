"""
Tests for EPyR Tools configuration system.
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from epyr.config import (
    EPyRConfig,
    config,
    get_performance_config,
    get_plotting_config,
    is_cache_enabled,
    is_debug_mode,
)


class TestEPyRConfig:
    """Test EPyRConfig class functionality."""

    def test_config_initialization(self):
        """Test config initialization with defaults."""
        test_config = EPyRConfig()

        # Check that defaults are loaded
        assert test_config.get("plotting.dpi") == 300
        assert test_config.get("performance.cache_enabled") is True
        assert test_config.get("logging.level") == "INFO"

    def test_get_nested_values(self):
        """Test getting nested configuration values."""
        test_config = EPyRConfig()

        # Test nested access
        assert test_config.get("plotting.dpi") == 300
        assert test_config.get("plotting.figure_size") == [8, 6]

        # Test non-existent keys
        assert test_config.get("nonexistent.key") is None
        assert test_config.get("nonexistent.key", "default") == "default"

    def test_set_nested_values(self):
        """Test setting nested configuration values."""
        test_config = EPyRConfig()

        # Set existing nested value
        test_config.set("plotting.dpi", 150)
        assert test_config.get("plotting.dpi") == 150

        # Set new nested value
        test_config.set("new.section.key", "value")
        assert test_config.get("new.section.key") == "value"

        # Set value that overwrites dict
        test_config.set("plotting", "not_a_dict")
        assert test_config.get("plotting") == "not_a_dict"

    def test_get_set_section(self):
        """Test getting and setting entire sections."""
        test_config = EPyRConfig()

        # Get existing section
        plotting_config = test_config.get_section("plotting")
        assert isinstance(plotting_config, dict)
        assert "dpi" in plotting_config

        # Set entire section
        new_section = {"key1": "value1", "key2": "value2"}
        test_config.set_section("new_section", new_section)
        assert test_config.get_section("new_section") == new_section

        # Get non-existent section
        empty_section = test_config.get_section("nonexistent")
        assert empty_section == {}

    def test_reset_section(self):
        """Test resetting section to defaults."""
        test_config = EPyRConfig()

        # Modify a value
        test_config.set("plotting.dpi", 999)
        assert test_config.get("plotting.dpi") == 999

        # Reset section
        test_config.reset_section("plotting")
        assert test_config.get("plotting.dpi") == 300  # Back to default

    def test_reset_all(self):
        """Test resetting all configuration to defaults."""
        test_config = EPyRConfig()

        # Modify several values
        test_config.set("plotting.dpi", 999)
        test_config.set("performance.cache_enabled", False)
        test_config.set("new.key", "value")

        # Reset all
        test_config.reset_all()

        # Check defaults are restored
        assert test_config.get("plotting.dpi") == 300
        assert test_config.get("performance.cache_enabled") is True
        assert test_config.get("new.key") is None

    def test_config_file_path(self):
        """Test configuration file path generation."""
        test_config = EPyRConfig()

        config_file = test_config.get_config_file_path()
        assert isinstance(config_file, Path)
        assert config_file.name == "config.json"

        # Should be in appropriate config directory
        config_dir = config_file.parent
        assert "epyrtools" in str(config_dir).lower()

    def test_save_load_config_file(self):
        """Test saving and loading configuration from file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            config_file = Path(temp_dir) / "test_config.json"

            # Create config with custom file path
            test_config = EPyRConfig()
            test_config._config_file = config_file

            # Modify some values
            test_config.set("plotting.dpi", 150)
            test_config.set("custom.key", "test_value")

            # Save to file
            test_config.save()
            assert config_file.exists()

            # Create new config instance and load from file
            new_config = EPyRConfig()
            new_config._config_file = config_file
            new_config._load_from_file()

            # Check values were loaded
            assert new_config.get("plotting.dpi") == 150
            assert new_config.get("custom.key") == "test_value"

    def test_export_import_config(self):
        """Test exporting and importing configuration."""
        with tempfile.TemporaryDirectory() as temp_dir:
            export_file = Path(temp_dir) / "exported_config.json"

            test_config = EPyRConfig()
            test_config.set("plotting.dpi", 150)
            test_config.set("export.test", "value")

            # Export configuration
            test_config.export_config(export_file)
            assert export_file.exists()

            # Check exported content
            with open(export_file, "r") as f:
                exported_data = json.load(f)
            assert exported_data["plotting"]["dpi"] == 150
            assert exported_data["export"]["test"] == "value"

            # Import into new config
            new_config = EPyRConfig()
            new_config.import_config(export_file)

            assert new_config.get("plotting.dpi") == 150
            assert new_config.get("export.test") == "value"

    def test_import_nonexistent_config(self):
        """Test importing from non-existent file."""
        test_config = EPyRConfig()

        with pytest.raises(FileNotFoundError):
            test_config.import_config(Path("/nonexistent/config.json"))

    def test_environment_variable_loading(self):
        """Test loading configuration from environment variables."""
        with patch.dict(
            os.environ,
            {
                "EPYR_PLOTTING_DPI": "150",
                "EPYR_PERFORMANCE_CACHE_ENABLED": "false",
                "EPYR_CUSTOM_KEY": '{"nested": "value"}',
            },
        ):
            test_config = EPyRConfig()

            # Check environment variables were loaded
            assert test_config.get("plotting.dpi") == "150"  # String from env
            assert (
                test_config.get("performance.cache_enabled") == "false"
            )  # String from env
            assert test_config.get("custom.key") == {"nested": "value"}  # JSON parsed

    def test_config_file_error_handling(self):
        """Test error handling for config file operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create invalid JSON file
            invalid_config_file = Path(temp_dir) / "invalid.json"
            invalid_config_file.write_text("invalid json content")

            test_config = EPyRConfig()
            test_config._config_file = invalid_config_file

            # Should not raise exception, just log warning
            test_config._load_from_file()

            # Should still have defaults
            assert test_config.get("plotting.dpi") == 300

    def test_repr(self):
        """Test string representation of config."""
        test_config = EPyRConfig()
        repr_str = repr(test_config)

        assert "EPyRConfig" in repr_str
        assert "config.json" in repr_str


class TestGlobalConfig:
    """Test global configuration instance."""

    def test_global_config_instance(self):
        """Test that config is properly initialized."""
        assert config is not None
        assert isinstance(config, EPyRConfig)

        # Test basic functionality
        assert config.get("plotting.dpi") == 300
        assert isinstance(config.get("plotting.figure_size"), list)

    def test_config_persistence_across_imports(self):
        """Test that config changes persist across imports."""
        # Modify config
        original_dpi = config.get("plotting.dpi")
        config.set("plotting.dpi", 999)

        try:
            # Re-import config
            from epyr.config import config as config2

            # Should be the same instance
            assert config is config2
            assert config2.get("plotting.dpi") == 999

        finally:
            # Restore original value
            config.set("plotting.dpi", original_dpi)


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_get_plotting_config(self):
        """Test get_plotting_config function."""
        plotting_config = get_plotting_config()

        assert isinstance(plotting_config, dict)
        assert "dpi" in plotting_config
        assert "default_style" in plotting_config
        assert "figure_size" in plotting_config

    def test_get_performance_config(self):
        """Test get_performance_config function."""
        performance_config = get_performance_config()

        assert isinstance(performance_config, dict)
        assert "cache_enabled" in performance_config
        assert "cache_size_mb" in performance_config
        assert "memory_limit_mb" in performance_config

    def test_is_debug_mode(self):
        """Test is_debug_mode function."""
        original_value = config.get("advanced.debug_mode")

        try:
            # Test default (should be False)
            assert is_debug_mode() is False

            # Set to True
            config.set("advanced.debug_mode", True)
            assert is_debug_mode() is True

        finally:
            config.set("advanced.debug_mode", original_value)

    def test_is_cache_enabled(self):
        """Test is_cache_enabled function."""
        original_value = config.get("performance.cache_enabled")

        try:
            # Test default (should be True)
            assert is_cache_enabled() is True

            # Set to False
            config.set("performance.cache_enabled", False)
            assert is_cache_enabled() is False

        finally:
            config.set("performance.cache_enabled", original_value)


class TestConfigDefaults:
    """Test default configuration values."""

    def test_plotting_defaults(self):
        """Test plotting section defaults."""
        test_config = EPyRConfig()

        plotting = test_config.get_section("plotting")

        assert plotting["default_style"] == "publication"
        assert plotting["dpi"] == 300
        assert plotting["figure_size"] == [8, 6]
        assert plotting["color_scheme"] == "viridis"
        assert plotting["font_size"] == 12
        assert plotting["line_width"] == 1.5
        assert plotting["grid_alpha"] == 0.3
        assert plotting["save_format"] == "png"

    def test_data_loading_defaults(self):
        """Test data loading section defaults."""
        test_config = EPyRConfig()

        data_loading = test_config.get_section("data_loading")

        assert data_loading["auto_plot"] is True
        assert data_loading["scaling_default"] == ""
        assert data_loading["file_dialog_remember_dir"] is True
        assert ".dta" in data_loading["supported_extensions"]
        assert ".dsc" in data_loading["supported_extensions"]

    def test_baseline_defaults(self):
        """Test baseline section defaults."""
        test_config = EPyRConfig()

        baseline = test_config.get_section("baseline")

        assert baseline["default_poly_order"] == 1
        assert baseline["default_method"] == "polynomial"
        assert baseline["exclusion_buffer"] == 0.1
        assert baseline["max_iterations"] == 1000

    def test_fair_conversion_defaults(self):
        """Test FAIR conversion section defaults."""
        test_config = EPyRConfig()

        fair = test_config.get_section("fair_conversion")

        assert fair["default_formats"] == ["csv", "json"]
        assert fair["include_metadata"] is True
        assert fair["preserve_precision"] is True
        assert fair["compression"] == "gzip"

    def test_performance_defaults(self):
        """Test performance section defaults."""
        test_config = EPyRConfig()

        performance = test_config.get_section("performance")

        assert performance["cache_enabled"] is True
        assert performance["cache_size_mb"] == 100
        assert performance["chunk_size_mb"] == 10
        assert performance["memory_limit_mb"] == 500
        assert performance["parallel_processing"] is True

    def test_logging_defaults(self):
        """Test logging section defaults."""
        test_config = EPyRConfig()

        logging_config = test_config.get_section("logging")

        assert logging_config["level"] == "INFO"
        assert logging_config["file_logging"] is False
        assert logging_config["log_file"] is None
        assert logging_config["console_output"] is True

    def test_gui_defaults(self):
        """Test GUI section defaults."""
        test_config = EPyRConfig()

        gui = test_config.get_section("gui")

        assert gui["theme"] == "default"
        assert gui["window_size"] == [800, 600]
        assert gui["remember_position"] is True
        assert gui["auto_refresh"] is True

    def test_advanced_defaults(self):
        """Test advanced section defaults."""
        test_config = EPyRConfig()

        advanced = test_config.get_section("advanced")

        assert advanced["debug_mode"] is False
        assert advanced["developer_mode"] is False
        assert advanced["experimental_features"] is False
        assert advanced["error_reporting"] is True


class TestConfigValidation:
    """Test configuration validation and error handling."""

    def test_deep_nested_access(self):
        """Test deeply nested configuration access."""
        test_config = EPyRConfig()

        # Set deeply nested value
        test_config.set("level1.level2.level3.key", "deep_value")
        assert test_config.get("level1.level2.level3.key") == "deep_value"

        # Access partial paths
        assert isinstance(test_config.get("level1"), dict)
        assert isinstance(test_config.get("level1.level2"), dict)
        assert isinstance(test_config.get("level1.level2.level3"), dict)

    def test_config_type_handling(self):
        """Test handling of different data types in configuration."""
        test_config = EPyRConfig()

        # Test different types
        test_config.set("test.string", "hello")
        test_config.set("test.integer", 42)
        test_config.set("test.float", 3.14)
        test_config.set("test.boolean", True)
        test_config.set("test.list", [1, 2, 3])
        test_config.set("test.dict", {"nested": "value"})

        assert test_config.get("test.string") == "hello"
        assert test_config.get("test.integer") == 42
        assert test_config.get("test.float") == 3.14
        assert test_config.get("test.boolean") is True
        assert test_config.get("test.list") == [1, 2, 3]
        assert test_config.get("test.dict") == {"nested": "value"}

    def test_config_merge_behavior(self):
        """Test configuration merging behavior."""
        test_config = EPyRConfig()

        # Set initial nested structure
        test_config.set("merge.section", {"key1": "value1", "key2": "value2"})

        # Merge additional data
        additional_config = {
            "merge": {
                "section": {
                    "key2": "updated_value2",  # Should overwrite
                    "key3": "value3",  # Should add
                }
            }
        }

        test_config._merge_config(additional_config)

        # Check merge results
        assert test_config.get("merge.section.key1") == "value1"  # Preserved
        assert test_config.get("merge.section.key2") == "updated_value2"  # Updated
        assert test_config.get("merge.section.key3") == "value3"  # Added
