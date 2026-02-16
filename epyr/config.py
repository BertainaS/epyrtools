"""
Configuration Management for EPyR Tools
=======================================

Centralized configuration system for EPyR Tools with support for:
- Default settings
- User preferences
- Environment variables
- Configuration files

Usage:
    from epyr.config import config

    # Get configuration values
    plot_style = config.get('plotting.default_style')
    cache_enabled = config.get('performance.cache_enabled')

    # Set configuration values
    config.set('plotting.dpi', 300)
    config.save()
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Union


class EPyRConfig:
    """Centralized configuration manager for EPyR Tools."""

    # Default configuration values
    _defaults = {
        "plotting": {
            "default_style": "publication",
            "dpi": 300,
            "figure_size": [8, 6],
            "color_scheme": "viridis",
            "font_size": 12,
            "line_width": 1.5,
            "grid_alpha": 0.3,
            "save_format": "png",
        },
        "data_loading": {
            "auto_plot": True,
            "scaling_default": "",
            "file_dialog_remember_dir": True,
            "supported_extensions": [".dta", ".dsc", ".spc", ".par"],
        },
        "baseline": {
            "default_poly_order": 1,
            "default_method": "polynomial",
            "exclusion_buffer": 0.1,  # 10% buffer around peaks
            "max_iterations": 1000,
        },
        "fair_conversion": {
            "default_formats": ["csv", "json"],
            "include_metadata": True,
            "preserve_precision": True,
            "compression": "gzip",
        },
        "performance": {
            "cache_enabled": True,
            "cache_size_mb": 100,
            "chunk_size_mb": 10,
            "memory_limit_mb": 500,
            "parallel_processing": True,
        },
        "logging": {
            "level": "INFO",
            "file_logging": False,
            "log_file": None,
            "console_output": True,
        },
        "gui": {
            "theme": "default",
            "window_size": [800, 600],
            "remember_position": True,
            "auto_refresh": True,
        },
        "advanced": {
            "debug_mode": False,
            "developer_mode": False,
            "experimental_features": False,
            "error_reporting": True,
        },
    }

    def __init__(self):
        """Initialize configuration manager."""
        self._config = self._load_defaults()
        self._config_file = self._get_config_file_path()
        self._load_from_file()
        self._load_from_environment()

    def _load_defaults(self) -> Dict[str, Any]:
        """Load default configuration values."""
        import copy

        return copy.deepcopy(self._defaults)

    def _get_config_file_path(self) -> Path:
        """Get path to user configuration file."""
        # Try to find appropriate config directory
        config_dir = None

        if os.name == "nt":  # Windows
            config_dir = Path(os.environ.get("APPDATA", Path.home())) / "EPyRTools"
        else:  # Unix-like
            config_dir = (
                Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
                / "epyrtools"
            )

        config_dir.mkdir(parents=True, exist_ok=True)
        return config_dir / "config.json"

    def _load_from_file(self) -> None:
        """Load configuration from file if it exists."""
        if self._config_file.exists():
            try:
                with open(self._config_file, "r", encoding="utf-8") as f:
                    file_config = json.load(f)
                self._merge_config(file_config)
            except (json.JSONDecodeError, IOError) as e:
                from .logging_config import get_logger

                logger = get_logger(__name__)
                logger.warning(f"Failed to load config from {self._config_file}: {e}")

    def _load_from_environment(self) -> None:
        """Load configuration from environment variables."""
        # Environment variables with EPYR_ prefix override config
        for key, value in os.environ.items():
            if key.startswith("EPYR_"):
                config_key = key[5:].lower().replace("_", ".")
                try:
                    # Try to parse as JSON first (for complex values)
                    parsed_value = json.loads(value)
                except json.JSONDecodeError:
                    # Fallback to string value
                    parsed_value = value

                self._set_nested(config_key, parsed_value)

    def _merge_config(self, new_config: Dict[str, Any]) -> None:
        """Merge new configuration into existing config."""

        def merge_dicts(base: Dict, overlay: Dict) -> Dict:
            for key, value in overlay.items():
                if (
                    key in base
                    and isinstance(base[key], dict)
                    and isinstance(value, dict)
                ):
                    merge_dicts(base[key], value)
                else:
                    base[key] = value
            return base

        merge_dicts(self._config, new_config)

    def _set_nested(self, key_path: str, value: Any) -> None:
        """Set a nested configuration value using dot notation."""
        keys = key_path.split(".")
        current = self._config

        # Navigate to the parent dictionary
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            elif not isinstance(current[key], dict):
                # Convert to dict if not already
                current[key] = {}
            current = current[key]

        # Set the final value
        current[keys[-1]] = value

    def _get_nested(self, key_path: str) -> Any:
        """Get a nested configuration value using dot notation."""
        keys = key_path.split(".")
        current = self._config

        try:
            for key in keys:
                current = current[key]
            return current
        except (KeyError, TypeError):
            return None

    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value by key.

        Args:
            key: Configuration key (supports dot notation, e.g., 'plotting.dpi')
            default: Default value if key not found

        Returns:
            Configuration value or default
        """
        value = self._get_nested(key)
        return value if value is not None else default

    def set(self, key: str, value: Any) -> None:
        """Set configuration value by key.

        Args:
            key: Configuration key (supports dot notation)
            value: Value to set
        """
        self._set_nested(key, value)

    def get_section(self, section: str) -> Dict[str, Any]:
        """Get entire configuration section.

        Args:
            section: Section name (e.g., 'plotting')

        Returns:
            Dictionary of section configuration
        """
        return self.get(section, {})

    def set_section(self, section: str, config_dict: Dict[str, Any]) -> None:
        """Set entire configuration section.

        Args:
            section: Section name
            config_dict: Dictionary of configuration values
        """
        self.set(section, config_dict)

    def reset_section(self, section: str) -> None:
        """Reset section to default values.

        Args:
            section: Section name to reset
        """
        if section in self._defaults:
            self.set_section(section, self._defaults[section].copy())

    def reset_all(self) -> None:
        """Reset all configuration to defaults."""
        self._config = self._load_defaults()

    def save(self) -> None:
        """Save current configuration to file."""
        try:
            with open(self._config_file, "w", encoding="utf-8") as f:
                json.dump(self._config, f, indent=2, ensure_ascii=False)
        except IOError as e:
            from .logging_config import get_logger

            logger = get_logger(__name__)
            logger.error(f"Failed to save config to {self._config_file}: {e}")

    def get_config_file_path(self) -> Path:
        """Get path to configuration file."""
        return self._config_file

    def export_config(self, file_path: Union[str, Path]) -> None:
        """Export current configuration to specified file.

        Args:
            file_path: Path to export configuration
        """
        export_path = Path(file_path)
        try:
            with open(export_path, "w", encoding="utf-8") as f:
                json.dump(self._config, f, indent=2, ensure_ascii=False)
        except IOError as e:
            from .logging_config import get_logger

            logger = get_logger(__name__)
            logger.error(f"Failed to export config to {export_path}: {e}")

    def import_config(self, file_path: Union[str, Path]) -> None:
        """Import configuration from specified file.

        Args:
            file_path: Path to import configuration from
        """
        import_path = Path(file_path)
        if not import_path.exists():
            raise FileNotFoundError(f"Config file not found: {import_path}")

        try:
            with open(import_path, "r", encoding="utf-8") as f:
                imported_config = json.load(f)
            self._merge_config(imported_config)
        except (json.JSONDecodeError, IOError) as e:
            from .logging_config import get_logger

            logger = get_logger(__name__)
            logger.error(f"Failed to import config from {import_path}: {e}")
            raise

    def __repr__(self) -> str:
        """String representation of configuration."""
        return f"EPyRConfig(file={self._config_file})"


# Global configuration instance
config = EPyRConfig()


# Convenience functions for common operations
def get_plotting_config() -> Dict[str, Any]:
    """Get plotting configuration section."""
    return config.get_section("plotting")


def get_performance_config() -> Dict[str, Any]:
    """Get performance configuration section."""
    return config.get_section("performance")


def is_debug_mode() -> bool:
    """Check if debug mode is enabled."""
    return config.get("advanced.debug_mode", False)


def is_cache_enabled() -> bool:
    """Check if caching is enabled."""
    return config.get("performance.cache_enabled", True)
