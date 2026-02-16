"""
Plugin Architecture for EPyR Tools
==================================

Extensible plugin system for adding support for new file formats,
data processing methods, and export formats.

This module provides:
- Base plugin interfaces and abstract classes
- Plugin discovery and loading mechanism
- Format handler registration system
- Extension point management

Usage:
    # Register a new file format plugin
    from epyr.plugins import PluginManager, FileFormatPlugin

    class MyFormatPlugin(FileFormatPlugin):
        format_name = "myformat"
        file_extensions = [".myf", ".myformat"]

        def can_load(self, file_path: Path) -> bool:
            return file_path.suffix.lower() in self.file_extensions

        def load(self, file_path: Path) -> Tuple[np.ndarray, np.ndarray, dict]:
            # Implementation here
            pass

    # Register the plugin
    plugin_manager.register_plugin(MyFormatPlugin())
"""

import importlib
import inspect
import pkgutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Type, Union

import numpy as np

from .config import config
from .logging_config import get_logger

logger = get_logger(__name__)


class BasePlugin(ABC):
    """Base class for all EPyR Tools plugins."""

    # Plugin metadata
    plugin_name: str = ""
    plugin_version: str = "1.0.0"
    plugin_description: str = ""
    plugin_author: str = ""

    def __init__(self):
        """Initialize the plugin."""
        if not self.plugin_name:
            self.plugin_name = self.__class__.__name__

    @abstractmethod
    def initialize(self) -> bool:
        """Initialize the plugin. Called when plugin is loaded.

        Returns:
            True if initialization successful, False otherwise
        """
        pass

    def cleanup(self):
        """Cleanup plugin resources. Called when plugin is unloaded."""
        pass

    def get_info(self) -> Dict[str, Any]:
        """Get plugin information."""
        return {
            "name": self.plugin_name,
            "version": self.plugin_version,
            "description": self.plugin_description,
            "author": self.plugin_author,
            "class": self.__class__.__name__,
        }


class FileFormatPlugin(BasePlugin):
    """Base class for file format plugins."""

    format_name: str = ""
    file_extensions: List[str] = []
    supports_loading: bool = True
    supports_saving: bool = False

    @abstractmethod
    def can_load(self, file_path: Path) -> bool:
        """Check if this plugin can load the given file.

        Args:
            file_path: Path to the file to check

        Returns:
            True if plugin can load this file
        """
        pass

    @abstractmethod
    def load(self, file_path: Path, **kwargs) -> Tuple[
        Optional[np.ndarray],
        Optional[Union[np.ndarray, List[np.ndarray]]],
        Optional[Dict[str, Any]],
    ]:
        """Load data from file.

        Args:
            file_path: Path to file to load
            **kwargs: Additional loading options

        Returns:
            Tuple of (x_data, y_data, parameters)
        """
        pass

    def can_save(self, data_format: str) -> bool:
        """Check if this plugin can save data in the specified format.

        Args:
            data_format: Format identifier

        Returns:
            True if plugin supports saving in this format
        """
        return self.supports_saving and data_format.lower() == self.format_name.lower()

    def save(
        self,
        file_path: Path,
        x_data: np.ndarray,
        y_data: np.ndarray,
        parameters: Dict[str, Any],
        **kwargs,
    ) -> bool:
        """Save data to file.

        Args:
            file_path: Path where to save the file
            x_data: X-axis data
            y_data: Y-axis data
            parameters: Metadata parameters
            **kwargs: Additional saving options

        Returns:
            True if save successful
        """
        if not self.supports_saving:
            raise NotImplementedError(
                f"Plugin {self.plugin_name} does not support saving"
            )
        return False


class ProcessingPlugin(BasePlugin):
    """Base class for data processing plugins."""

    processing_name: str = ""
    input_requirements: List[str] = []  # e.g., ['1d_data', '2d_data']
    output_types: List[str] = []  # e.g., ['corrected_data', 'fit_parameters']

    @abstractmethod
    def process(self, data: Dict[str, Any], **kwargs) -> Dict[str, Any]:
        """Process EPR data.

        Args:
            data: Input data dictionary
            **kwargs: Processing parameters

        Returns:
            Dictionary with processed results
        """
        pass

    def validate_input(self, data: Dict[str, Any]) -> bool:
        """Validate input data requirements.

        Args:
            data: Input data dictionary

        Returns:
            True if input data meets requirements
        """
        for requirement in self.input_requirements:
            if requirement not in data:
                return False
        return True


class ExportPlugin(BasePlugin):
    """Base class for data export plugins."""

    export_format: str = ""
    file_extension: str = ""
    supports_metadata: bool = True

    @abstractmethod
    def export(
        self,
        output_path: Path,
        x_data: np.ndarray,
        y_data: np.ndarray,
        parameters: Dict[str, Any],
        **kwargs,
    ) -> bool:
        """Export data to specified format.

        Args:
            output_path: Path for output file
            x_data: X-axis data
            y_data: Y-axis data
            parameters: Metadata parameters
            **kwargs: Export options

        Returns:
            True if export successful
        """
        pass


class PluginManager:
    """Manages plugin discovery, loading, and registration."""

    def __init__(self):
        """Initialize plugin manager."""
        self.file_format_plugins: Dict[str, FileFormatPlugin] = {}
        self.processing_plugins: Dict[str, ProcessingPlugin] = {}
        self.export_plugins: Dict[str, ExportPlugin] = {}
        self.loaded_plugins: Dict[str, BasePlugin] = {}

    def register_plugin(self, plugin: BasePlugin) -> bool:
        """Register a plugin instance.

        Args:
            plugin: Plugin instance to register

        Returns:
            True if registration successful
        """
        try:
            # Initialize the plugin
            if not plugin.initialize():
                logger.error(f"Failed to initialize plugin {plugin.plugin_name}")
                return False

            # Register based on plugin type
            if isinstance(plugin, FileFormatPlugin):
                self.file_format_plugins[plugin.format_name] = plugin
                logger.info(f"Registered file format plugin: {plugin.format_name}")
            elif isinstance(plugin, ProcessingPlugin):
                self.processing_plugins[plugin.processing_name] = plugin
                logger.info(f"Registered processing plugin: {plugin.processing_name}")
            elif isinstance(plugin, ExportPlugin):
                self.export_plugins[plugin.export_format] = plugin
                logger.info(f"Registered export plugin: {plugin.export_format}")

            self.loaded_plugins[plugin.plugin_name] = plugin
            return True

        except Exception as e:
            logger.error(f"Failed to register plugin {plugin.plugin_name}: {e}")
            return False

    def unregister_plugin(self, plugin_name: str) -> bool:
        """Unregister a plugin by name.

        Args:
            plugin_name: Name of plugin to unregister

        Returns:
            True if unregistration successful
        """
        if plugin_name not in self.loaded_plugins:
            logger.warning(f"Plugin {plugin_name} not found")
            return False

        plugin = self.loaded_plugins[plugin_name]

        try:
            # Cleanup plugin
            plugin.cleanup()

            # Remove from specific registries
            if isinstance(plugin, FileFormatPlugin):
                self.file_format_plugins.pop(plugin.format_name, None)
            elif isinstance(plugin, ProcessingPlugin):
                self.processing_plugins.pop(plugin.processing_name, None)
            elif isinstance(plugin, ExportPlugin):
                self.export_plugins.pop(plugin.export_format, None)

            # Remove from loaded plugins
            del self.loaded_plugins[plugin_name]
            logger.info(f"Unregistered plugin: {plugin_name}")
            return True

        except Exception as e:
            logger.error(f"Failed to unregister plugin {plugin_name}: {e}")
            return False

    def discover_plugins(self, plugin_directories: Optional[List[Path]] = None) -> int:
        """Discover and load plugins from specified directories.

        Args:
            plugin_directories: List of directories to search for plugins

        Returns:
            Number of plugins loaded
        """
        if plugin_directories is None:
            plugin_directories = self._get_default_plugin_directories()

        loaded_count = 0

        for plugin_dir in plugin_directories:
            if not plugin_dir.exists():
                logger.debug(f"Plugin directory does not exist: {plugin_dir}")
                continue

            logger.debug(f"Scanning for plugins in: {plugin_dir}")

            try:
                # Add directory to Python path temporarily
                import sys

                if str(plugin_dir) not in sys.path:
                    sys.path.insert(0, str(plugin_dir))

                # Discover Python modules in directory
                for finder, name, ispkg in pkgutil.iter_modules([str(plugin_dir)]):
                    try:
                        module = importlib.import_module(name)
                        plugins_found = self._extract_plugins_from_module(module)
                        loaded_count += plugins_found

                    except Exception as e:
                        logger.warning(f"Failed to load plugin module {name}: {e}")
                        continue

            except Exception as e:
                logger.error(f"Error scanning plugin directory {plugin_dir}: {e}")
                continue

        logger.info(f"Plugin discovery complete. Loaded {loaded_count} plugins.")
        return loaded_count

    def _get_default_plugin_directories(self) -> List[Path]:
        """Get default directories to search for plugins."""
        directories = []

        # User plugin directory
        config_dir = config.get_config_file_path().parent
        user_plugin_dir = config_dir / "plugins"
        directories.append(user_plugin_dir)

        # System plugin directory
        try:
            import epyr

            package_dir = Path(epyr.__file__).parent
            system_plugin_dir = package_dir / "plugins"
            directories.append(system_plugin_dir)
        except:
            pass

        return directories

    def _extract_plugins_from_module(self, module) -> int:
        """Extract plugin classes from a module and register them.

        Args:
            module: Python module to examine

        Returns:
            Number of plugins found and registered
        """
        loaded_count = 0

        for name, obj in inspect.getmembers(module):
            if (
                inspect.isclass(obj)
                and issubclass(obj, BasePlugin)
                and obj is not BasePlugin
                and obj not in [FileFormatPlugin, ProcessingPlugin, ExportPlugin]
            ):

                try:
                    plugin_instance = obj()
                    if self.register_plugin(plugin_instance):
                        loaded_count += 1
                except Exception as e:
                    logger.warning(f"Failed to instantiate plugin {name}: {e}")

        return loaded_count

    def get_file_format_plugin(self, file_path: Path) -> Optional[FileFormatPlugin]:
        """Get appropriate file format plugin for a file.

        Args:
            file_path: Path to file

        Returns:
            Plugin that can load the file, or None
        """
        for plugin in self.file_format_plugins.values():
            if plugin.can_load(file_path):
                return plugin
        return None

    def get_export_plugin(self, format_name: str) -> Optional[ExportPlugin]:
        """Get export plugin for specified format.

        Args:
            format_name: Export format name

        Returns:
            Plugin that can export to format, or None
        """
        return self.export_plugins.get(format_name.lower())

    def get_processing_plugin(self, processing_name: str) -> Optional[ProcessingPlugin]:
        """Get processing plugin by name.

        Args:
            processing_name: Processing method name

        Returns:
            Processing plugin, or None
        """
        return self.processing_plugins.get(processing_name.lower())

    def list_plugins(self) -> Dict[str, List[Dict[str, Any]]]:
        """List all loaded plugins by category.

        Returns:
            Dictionary with plugin information by category
        """
        return {
            "file_formats": [
                plugin.get_info() for plugin in self.file_format_plugins.values()
            ],
            "processing": [
                plugin.get_info() for plugin in self.processing_plugins.values()
            ],
            "export": [plugin.get_info() for plugin in self.export_plugins.values()],
        }

    def get_supported_extensions(self) -> List[str]:
        """Get list of all supported file extensions.

        Returns:
            List of file extensions (including dots)
        """
        extensions = []
        for plugin in self.file_format_plugins.values():
            extensions.extend(plugin.file_extensions)
        return sorted(list(set(extensions)))


# Example built-in plugins


class CSVExportPlugin(ExportPlugin):
    """Built-in CSV export plugin."""

    plugin_name = "CSV Exporter"
    export_format = "csv"
    file_extension = ".csv"

    def initialize(self) -> bool:
        return True

    def export(
        self,
        output_path: Path,
        x_data: np.ndarray,
        y_data: np.ndarray,
        parameters: Dict[str, Any],
        **kwargs,
    ) -> bool:
        """Export data to CSV format."""
        try:
            import pandas as pd

            # Create DataFrame
            df = pd.DataFrame(
                {
                    "field": (
                        x_data if hasattr(x_data, "__len__") else range(len(y_data))
                    ),
                    "intensity": y_data,
                }
            )

            # Add metadata as comments if supported
            metadata_comment = f"# Parameters: {parameters}" if parameters else ""

            # Save to CSV
            df.to_csv(output_path, index=False)

            # Add metadata comment at the top if parameters provided
            if metadata_comment and kwargs.get("include_metadata", True):
                content = output_path.read_text()
                output_path.write_text(f"{metadata_comment}\n{content}")

            return True

        except Exception as e:
            logger.error(f"CSV export failed: {e}")
            return False


# Global plugin manager instance
plugin_manager = PluginManager()

# Register built-in plugins
plugin_manager.register_plugin(CSVExportPlugin())

# Auto-discover plugins on import
if config.get("advanced.experimental_features", False):
    plugin_manager.discover_plugins()
