"""
Tests for EPyR Tools plugin system.
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

from epyr.plugins import (
    BasePlugin,
    ExportPlugin,
    FileFormatPlugin,
    PluginManager,
    ProcessingPlugin,
)


class MockFileFormatPlugin(FileFormatPlugin):
    """Mock file format plugin for testing."""

    plugin_name = "Mock Format Plugin"
    format_name = "mockfmt"
    file_extensions = [".mock", ".test"]
    supports_loading = True

    def __init__(self):
        super().__init__()
        self.initialized = False
        self.load_called = False

    def initialize(self) -> bool:
        self.initialized = True
        return True

    def can_load(self, file_path: Path) -> bool:
        return file_path.suffix.lower() in self.file_extensions

    def load(self, file_path: Path, **kwargs):
        self.load_called = True
        # Return mock EPR data
        x_data = np.linspace(3400, 3500, 100)
        y_data = np.sin(x_data)
        params = {"frequency": 9.4e9, "file_path": str(file_path)}
        return x_data, y_data, params


class MockProcessingPlugin(ProcessingPlugin):
    """Mock processing plugin for testing."""

    plugin_name = "Mock Processing Plugin"
    processing_name = "mockprocess"
    input_requirements = ["y_data"]
    output_types = ["processed_data"]

    def __init__(self):
        super().__init__()
        self.initialized = False
        self.process_called = False

    def initialize(self) -> bool:
        self.initialized = True
        return True

    def process(self, data: dict, **kwargs) -> dict:
        self.process_called = True
        # Simple processing: multiply by 2
        result = data.copy()
        if "y_data" in data:
            result["processed_data"] = data["y_data"] * 2
        return result


class MockExportPlugin(ExportPlugin):
    """Mock export plugin for testing."""

    plugin_name = "Mock Export Plugin"
    export_format = "mockexport"
    file_extension = ".mockexp"
    supports_metadata = True

    def __init__(self):
        super().__init__()
        self.initialized = False
        self.export_called = False
        self.exported_files = []

    def initialize(self) -> bool:
        self.initialized = True
        return True

    def export(
        self, output_path: Path, x_data, y_data, parameters: dict, **kwargs
    ) -> bool:
        self.export_called = True
        self.exported_files.append(output_path)

        # Write mock export data
        with open(output_path, "w") as f:
            f.write(f"# Mock export\n")
            f.write(f"# Parameters: {parameters}\n")
            for x, y in zip(x_data[:5], y_data[:5]):  # First 5 points
                f.write(f"{x},{y}\n")

        return True


class FailingPlugin(BasePlugin):
    """Plugin that fails during initialization."""

    plugin_name = "Failing Plugin"

    def initialize(self) -> bool:
        raise Exception("Initialization failed")


class TestBasePlugin:
    """Test base plugin functionality."""

    def test_base_plugin_creation(self):
        """Test base plugin instantiation."""

        class TestPlugin(BasePlugin):
            def initialize(self) -> bool:
                return True

        plugin = TestPlugin()
        assert plugin.plugin_name == "TestPlugin"
        assert plugin.plugin_version == "1.0.0"

        info = plugin.get_info()
        assert info["name"] == "TestPlugin"
        assert info["version"] == "1.0.0"
        assert info["class"] == "TestPlugin"

    def test_base_plugin_custom_name(self):
        """Test base plugin with custom name."""

        class TestPlugin(BasePlugin):
            plugin_name = "Custom Plugin Name"

            def initialize(self) -> bool:
                return True

        plugin = TestPlugin()
        assert plugin.plugin_name == "Custom Plugin Name"


class TestFileFormatPlugin:
    """Test file format plugin functionality."""

    def test_file_format_plugin_creation(self):
        """Test file format plugin instantiation."""
        plugin = MockFileFormatPlugin()

        assert plugin.format_name == "mockfmt"
        assert plugin.file_extensions == [".mock", ".test"]
        assert plugin.supports_loading is True
        assert plugin.supports_saving is False

    def test_can_load(self):
        """Test file format detection."""
        plugin = MockFileFormatPlugin()

        assert plugin.can_load(Path("test.mock")) is True
        assert plugin.can_load(Path("test.TEST")) is True  # Case insensitive
        assert plugin.can_load(Path("test.dsc")) is False

    def test_load_functionality(self):
        """Test data loading."""
        plugin = MockFileFormatPlugin()

        with tempfile.NamedTemporaryFile(suffix=".mock", delete=False) as f:
            test_file = Path(f.name)

        try:
            x_data, y_data, params = plugin.load(test_file)

            assert plugin.load_called is True
            assert isinstance(x_data, np.ndarray)
            assert isinstance(y_data, np.ndarray)
            assert isinstance(params, dict)
            assert "frequency" in params
        finally:
            test_file.unlink()

    def test_save_not_supported(self):
        """Test that saving raises error when not supported."""
        plugin = MockFileFormatPlugin()

        with pytest.raises(NotImplementedError):
            plugin.save(Path("test.mock"), None, None, {})

    def test_can_save(self):
        """Test save capability detection."""
        plugin = MockFileFormatPlugin()

        # Default plugin doesn't support saving
        assert plugin.can_save("mockfmt") is False
        assert plugin.can_save("csv") is False


class TestProcessingPlugin:
    """Test processing plugin functionality."""

    def test_processing_plugin_creation(self):
        """Test processing plugin instantiation."""
        plugin = MockProcessingPlugin()

        assert plugin.processing_name == "mockprocess"
        assert plugin.input_requirements == ["y_data"]
        assert plugin.output_types == ["processed_data"]

    def test_validate_input_success(self):
        """Test input validation with valid data."""
        plugin = MockProcessingPlugin()

        data = {"y_data": np.array([1, 2, 3])}
        assert plugin.validate_input(data) is True

    def test_validate_input_failure(self):
        """Test input validation with invalid data."""
        plugin = MockProcessingPlugin()

        data = {"x_data": np.array([1, 2, 3])}  # Missing y_data
        assert plugin.validate_input(data) is False

    def test_process_functionality(self):
        """Test data processing."""
        plugin = MockProcessingPlugin()

        input_data = {"y_data": np.array([1, 2, 3])}
        result = plugin.process(input_data)

        assert plugin.process_called is True
        assert "processed_data" in result
        np.testing.assert_array_equal(result["processed_data"], np.array([2, 4, 6]))


class TestExportPlugin:
    """Test export plugin functionality."""

    def test_export_plugin_creation(self):
        """Test export plugin instantiation."""
        plugin = MockExportPlugin()

        assert plugin.export_format == "mockexport"
        assert plugin.file_extension == ".mockexp"
        assert plugin.supports_metadata is True

    def test_export_functionality(self):
        """Test data export."""
        plugin = MockExportPlugin()

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_export.mockexp"
            x_data = np.linspace(1, 5, 5)
            y_data = np.array([1, 4, 9, 16, 25])
            parameters = {"test_param": "test_value"}

            result = plugin.export(output_path, x_data, y_data, parameters)

            assert result is True
            assert plugin.export_called is True
            assert output_path in plugin.exported_files
            assert output_path.exists()

            # Check file contents
            content = output_path.read_text()
            assert "Mock export" in content
            assert "test_value" in content


class TestPluginManager:
    """Test plugin manager functionality."""

    def test_plugin_manager_creation(self):
        """Test plugin manager instantiation."""
        manager = PluginManager()

        assert isinstance(manager.file_format_plugins, dict)
        assert isinstance(manager.processing_plugins, dict)
        assert isinstance(manager.export_plugins, dict)
        assert isinstance(manager.loaded_plugins, dict)

    def test_register_file_format_plugin(self):
        """Test registering file format plugin."""
        manager = PluginManager()
        plugin = MockFileFormatPlugin()

        success = manager.register_plugin(plugin)

        assert success is True
        assert plugin.initialized is True
        assert "mockfmt" in manager.file_format_plugins
        assert manager.file_format_plugins["mockfmt"] is plugin
        assert plugin.plugin_name in manager.loaded_plugins

    def test_register_processing_plugin(self):
        """Test registering processing plugin."""
        manager = PluginManager()
        plugin = MockProcessingPlugin()

        success = manager.register_plugin(plugin)

        assert success is True
        assert plugin.initialized is True
        assert "mockprocess" in manager.processing_plugins
        assert manager.processing_plugins["mockprocess"] is plugin

    def test_register_export_plugin(self):
        """Test registering export plugin."""
        manager = PluginManager()
        plugin = MockExportPlugin()

        success = manager.register_plugin(plugin)

        assert success is True
        assert plugin.initialized is True
        assert "mockexport" in manager.export_plugins
        assert manager.export_plugins["mockexport"] is plugin

    def test_register_failing_plugin(self):
        """Test registering plugin that fails initialization."""
        manager = PluginManager()
        plugin = FailingPlugin()

        success = manager.register_plugin(plugin)

        assert success is False
        assert plugin.plugin_name not in manager.loaded_plugins

    def test_unregister_plugin(self):
        """Test unregistering plugin."""
        manager = PluginManager()
        plugin = MockFileFormatPlugin()

        # Register first
        manager.register_plugin(plugin)
        assert plugin.plugin_name in manager.loaded_plugins

        # Unregister
        success = manager.unregister_plugin(plugin.plugin_name)

        assert success is True
        assert plugin.plugin_name not in manager.loaded_plugins
        assert "mockfmt" not in manager.file_format_plugins

    def test_unregister_nonexistent_plugin(self):
        """Test unregistering plugin that doesn't exist."""
        manager = PluginManager()

        success = manager.unregister_plugin("nonexistent")

        assert success is False

    def test_get_file_format_plugin(self):
        """Test finding file format plugin."""
        manager = PluginManager()
        plugin = MockFileFormatPlugin()
        manager.register_plugin(plugin)

        # Test finding plugin
        found_plugin = manager.get_file_format_plugin(Path("test.mock"))
        assert found_plugin is plugin

        # Test not finding plugin
        not_found = manager.get_file_format_plugin(Path("test.unknown"))
        assert not_found is None

    def test_get_export_plugin(self):
        """Test finding export plugin."""
        manager = PluginManager()
        plugin = MockExportPlugin()
        manager.register_plugin(plugin)

        # Test finding plugin
        found_plugin = manager.get_export_plugin("mockexport")
        assert found_plugin is plugin

        # Test case insensitive
        found_plugin = manager.get_export_plugin("MOCKEXPORT")
        assert found_plugin is plugin

        # Test not finding plugin
        not_found = manager.get_export_plugin("unknown")
        assert not_found is None

    def test_get_processing_plugin(self):
        """Test finding processing plugin."""
        manager = PluginManager()
        plugin = MockProcessingPlugin()
        manager.register_plugin(plugin)

        # Test finding plugin
        found_plugin = manager.get_processing_plugin("mockprocess")
        assert found_plugin is plugin

        # Test not finding plugin
        not_found = manager.get_processing_plugin("unknown")
        assert not_found is None

    def test_list_plugins(self):
        """Test listing all plugins."""
        manager = PluginManager()

        file_plugin = MockFileFormatPlugin()
        processing_plugin = MockProcessingPlugin()
        export_plugin = MockExportPlugin()

        manager.register_plugin(file_plugin)
        manager.register_plugin(processing_plugin)
        manager.register_plugin(export_plugin)

        plugins_list = manager.list_plugins()

        assert "file_formats" in plugins_list
        assert "processing" in plugins_list
        assert "export" in plugins_list

        assert len(plugins_list["file_formats"]) == 1
        assert len(plugins_list["processing"]) == 1
        assert len(plugins_list["export"]) == 1

        # Check plugin info structure
        file_info = plugins_list["file_formats"][0]
        assert file_info["name"] == "Mock Format Plugin"
        assert file_info["class"] == "MockFileFormatPlugin"

    def test_get_supported_extensions(self):
        """Test getting supported file extensions."""
        manager = PluginManager()
        plugin = MockFileFormatPlugin()
        manager.register_plugin(plugin)

        extensions = manager.get_supported_extensions()

        assert ".mock" in extensions
        assert ".test" in extensions
        assert isinstance(extensions, list)
        assert extensions == sorted(extensions)  # Should be sorted

    def test_discover_plugins_nonexistent_dir(self):
        """Test plugin discovery with nonexistent directory."""
        manager = PluginManager()

        count = manager.discover_plugins([Path("/nonexistent/directory")])

        assert count == 0

    @patch("pkgutil.iter_modules")
    @patch("importlib.import_module")
    def test_discover_plugins_mock(self, mock_import, mock_iter):
        """Test plugin discovery with mocked modules."""
        manager = PluginManager()

        # Mock module discovery
        mock_iter.return_value = [(None, "test_plugin", False)]

        # Mock module with plugin class
        mock_module = Mock()
        mock_module.TestPlugin = MockFileFormatPlugin
        mock_import.return_value = mock_module

        with tempfile.TemporaryDirectory() as temp_dir:
            count = manager.discover_plugins([Path(temp_dir)])

        assert count == 1
        assert "mockfmt" in manager.file_format_plugins


class TestBuiltinPlugins:
    """Test built-in plugins."""

    def test_csv_export_plugin_registered(self):
        """Test that CSV export plugin is registered by default."""
        # Import plugin manager to trigger auto-registration
        from epyr.plugins import plugin_manager

        csv_plugin = plugin_manager.get_export_plugin("csv")
        assert csv_plugin is not None
        assert csv_plugin.export_format == "csv"

    def test_csv_export_functionality(self):
        """Test CSV export plugin functionality."""
        from epyr.plugins import plugin_manager

        csv_plugin = plugin_manager.get_export_plugin("csv")

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test.csv"
            x_data = np.array([1, 2, 3, 4, 5])
            y_data = np.array([1, 4, 9, 16, 25])
            parameters = {"test": "value"}

            success = csv_plugin.export(output_path, x_data, y_data, parameters)

            assert success is True
            assert output_path.exists()

            # Check CSV content
            content = output_path.read_text()
            assert "field,intensity" in content  # Header
            assert "1,1" in content  # Data points


class TestPluginIntegration:
    """Test plugin integration with EPyR Tools."""

    def test_plugin_manager_singleton(self):
        """Test that plugin_manager is a singleton instance."""
        # Import again to ensure same instance
        from epyr.plugins import plugin_manager
        from epyr.plugins import plugin_manager as manager2

        assert plugin_manager is manager2

    def test_plugin_configuration_integration(self):
        """Test plugin configuration with EPyR config system."""
        from epyr.config import config
        from epyr.plugins import plugin_manager

        # Test experimental features flag
        original_value = config.get("advanced.experimental_features")

        try:
            # Enable experimental features
            config.set("advanced.experimental_features", True)

            # This would trigger plugin discovery in real scenario
            # Here we just verify the config is accessible
            assert config.get("advanced.experimental_features") is True

        finally:
            # Restore original value
            config.set("advanced.experimental_features", original_value)
