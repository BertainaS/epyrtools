"""
Tests for EPyR Tools CLI functionality.
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from epyr import cli
from epyr.config import config


class TestCLICommands:
    """Test CLI command functions."""

    def test_cmd_info_basic(self, capsys):
        """Test basic info command functionality."""
        with patch("sys.argv", ["epyr-info"]):
            cli.cmd_info()

        captured = capsys.readouterr()
        assert "EPyR Tools Version" in captured.out
        assert "Configuration file" in captured.out

    def test_cmd_info_config(self, capsys):
        """Test info command with config flag."""
        with patch("sys.argv", ["epyr-info", "--config"]):
            cli.cmd_info()

        captured = capsys.readouterr()
        assert "Configuration" in captured.out
        assert "plotting" in captured.out  # Should show config sections

    def test_cmd_info_performance(self, capsys):
        """Test info command with performance flag."""
        with patch("sys.argv", ["epyr-info", "--performance"]):
            cli.cmd_info()

        captured = capsys.readouterr()
        assert "Performance Information" in captured.out
        assert "memory" in captured.out

    def test_cmd_info_plugins(self, capsys):
        """Test info command with plugins flag."""
        with patch("sys.argv", ["epyr-info", "--plugins"]):
            cli.cmd_info()

        captured = capsys.readouterr()
        assert "Loaded Plugins" in captured.out

    def test_cmd_info_all(self, capsys):
        """Test info command with all flags."""
        with patch("sys.argv", ["epyr-info", "--all"]):
            cli.cmd_info()

        captured = capsys.readouterr()
        assert "Configuration" in captured.out
        assert "Performance Information" in captured.out
        assert "Loaded Plugins" in captured.out


class TestConfigCommand:
    """Test configuration CLI commands."""

    def test_cmd_config_show(self, capsys):
        """Test config show command."""
        with patch("sys.argv", ["epyr-config", "show"]):
            cli.cmd_config()

        captured = capsys.readouterr()
        assert "plotting" in captured.out
        assert "performance" in captured.out

    def test_cmd_config_show_section(self, capsys):
        """Test config show specific section."""
        with patch("sys.argv", ["epyr-config", "show", "plotting"]):
            cli.cmd_config()

        captured = capsys.readouterr()
        config_output = json.loads(captured.out)
        assert "dpi" in config_output
        assert "default_style" in config_output

    def test_cmd_config_set(self):
        """Test config set command."""
        original_value = config.get("plotting.dpi")

        try:
            with patch("sys.argv", ["epyr-config", "set", "plotting.dpi", "150"]):
                cli.cmd_config()

            # Check that value was set
            assert config.get("plotting.dpi") == 150
        finally:
            # Restore original value
            config.set("plotting.dpi", original_value)

    def test_cmd_config_set_json_value(self):
        """Test config set with JSON value."""
        original_value = config.get("plotting.figure_size")

        try:
            with patch(
                "sys.argv", ["epyr-config", "set", "plotting.figure_size", "[10, 8]"]
            ):
                cli.cmd_config()

            # Check that JSON value was parsed correctly
            assert config.get("plotting.figure_size") == [10, 8]
        finally:
            # Restore original value
            config.set("plotting.figure_size", original_value)

    def test_cmd_config_export_import(self):
        """Test config export and import."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_file = f.name

        try:
            # Export config
            with patch("sys.argv", ["epyr-config", "export", temp_file]):
                cli.cmd_config()

            # Check file was created
            assert Path(temp_file).exists()

            # Modify config
            original_dpi = config.get("plotting.dpi")
            config.set("plotting.dpi", 999)

            # Import config
            with patch("sys.argv", ["epyr-config", "import", temp_file]):
                cli.cmd_config()

            # Check config was restored
            assert config.get("plotting.dpi") == original_dpi

        finally:
            # Cleanup
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_cmd_config_reset_section(self):
        """Test config reset section."""
        # Modify a value
        original_value = config.get("plotting.dpi")
        config.set("plotting.dpi", 999)

        # Reset section
        with patch("sys.argv", ["epyr-config", "reset", "plotting"]):
            cli.cmd_config()

        # Check value was reset to default
        assert config.get("plotting.dpi") == 300  # Default value


class TestConvertCommand:
    """Test conversion CLI commands."""

    def test_cmd_convert_missing_file(self, capsys):
        """Test convert command with missing input file."""
        with patch("sys.argv", ["epyr-convert", "nonexistent.dsc"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.cmd_convert()

            assert excinfo.value.code == 1

    @patch("epyr.fair.convert_bruker_to_fair")
    def test_cmd_convert_success(self, mock_convert):
        """Test successful conversion."""
        mock_convert.return_value = True

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                with patch(
                    "sys.argv",
                    [
                        "epyr-convert",
                        test_file,
                        "--output-dir",
                        temp_dir,
                        "--formats",
                        "csv,json",
                    ],
                ):
                    cli.cmd_convert()

                # Check that conversion function was called
                mock_convert.assert_called_once()
                call_args = mock_convert.call_args
                assert call_args[1]["formats"] == ["csv", "json"]
                assert call_args[1]["output_dir"] == temp_dir
        finally:
            os.unlink(test_file)

    @patch("epyr.fair.convert_bruker_to_fair")
    def test_cmd_convert_failure(self, mock_convert, capsys):
        """Test conversion failure handling."""
        mock_convert.return_value = False

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with patch("sys.argv", ["epyr-convert", test_file]):
                with pytest.raises(SystemExit) as excinfo:
                    cli.cmd_convert()

                assert excinfo.value.code == 1
        finally:
            os.unlink(test_file)


class TestValidateCommand:
    """Test validation CLI commands."""

    def test_cmd_validate_missing_file(self, capsys):
        """Test validate command with missing files."""
        with patch("sys.argv", ["epyr-validate", "nonexistent.dsc"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.cmd_validate()

            assert excinfo.value.code == 1

    @patch("epyr.eprload")
    def test_cmd_validate_success(self, mock_eprload, capsys):
        """Test successful validation."""
        import numpy as np

        # Mock successful data loading
        mock_eprload.return_value = (
            np.linspace(3400, 3500, 100),  # x_data
            np.random.randn(100),  # y_data
            {"frequency": 9.4e9},  # params
            "test.dsc",  # file_path
        )

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with patch("sys.argv", ["epyr-validate", test_file]):
                cli.cmd_validate()

            captured = capsys.readouterr()
            assert "✓" in captured.out
            assert "Valid" in captured.out
            assert "1/1 files valid" in captured.out
        finally:
            os.unlink(test_file)

    @patch("epyr.eprload")
    def test_cmd_validate_detailed(self, mock_eprload, capsys):
        """Test detailed validation output."""
        import numpy as np

        # Mock successful data loading
        mock_eprload.return_value = (
            np.linspace(3400, 3500, 100),
            np.random.randn(100),
            {"frequency": 9.4e9},
            "test.dsc",
        )

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with patch("sys.argv", ["epyr-validate", test_file, "--detailed"]):
                cli.cmd_validate()

            captured = capsys.readouterr()
            assert "Data points:" in captured.out
            assert "X-axis range:" in captured.out
            assert "FAIR compliance:" in captured.out
        finally:
            os.unlink(test_file)

    @patch("epyr.eprload")
    def test_cmd_validate_failure(self, mock_eprload, capsys):
        """Test validation with invalid data."""
        # Mock failed data loading
        mock_eprload.return_value = (None, None, None, None)

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with patch("sys.argv", ["epyr-validate", test_file]):
                with pytest.raises(SystemExit) as excinfo:
                    cli.cmd_validate()

                assert excinfo.value.code == 1

            captured = capsys.readouterr()
            assert "✗" in captured.out
            assert "Invalid data" in captured.out
        finally:
            os.unlink(test_file)


class TestIsotopesCommand:
    """Test isotopes GUI command."""

    @patch("epyr.isotope_gui.run_gui")
    def test_cmd_isotopes_success(self, mock_run_gui):
        """Test successful isotopes GUI launch."""
        with patch("sys.argv", ["epyr-isotopes"]):
            cli.cmd_isotopes()

        mock_run_gui.assert_called_once()

    @patch("epyr.isotope_gui.run_gui")
    def test_cmd_isotopes_failure(self, mock_run_gui):
        """Test isotopes GUI launch failure."""
        mock_run_gui.side_effect = Exception("GUI failed to start")

        with patch("sys.argv", ["epyr-isotopes"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.cmd_isotopes()

            assert excinfo.value.code == 1


class TestMainCLI:
    """Test main CLI entry point."""

    def test_main_no_command(self, capsys):
        """Test main CLI with no command."""
        with patch("sys.argv", ["epyr"]):
            cli.main()

        captured = capsys.readouterr()
        assert "usage:" in captured.out
        assert "Available commands" in captured.out
        assert "convert" in captured.out

    def test_main_help(self, capsys):
        """Test main CLI help."""
        with patch("sys.argv", ["epyr", "--help"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.main()

            # Help should exit with code 0
            assert excinfo.value.code == 0

    @patch("epyr.cli.cmd_convert")
    def test_main_dispatch_convert(self, mock_cmd):
        """Test command dispatch to convert."""
        with patch("sys.argv", ["epyr", "convert"]):
            cli.main()

        mock_cmd.assert_called_once()

    @patch("epyr.cli.cmd_config")
    def test_main_dispatch_config(self, mock_cmd):
        """Test command dispatch to config."""
        with patch("sys.argv", ["epyr", "config"]):
            cli.main()

        mock_cmd.assert_called_once()


class TestBatchConvertCommand:
    """Test batch convert CLI command."""

    def test_cmd_batch_convert_missing_dir(self, capsys):
        """Test batch convert with missing directory."""
        with patch("sys.argv", ["epyr-batch-convert", "nonexistent_dir"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.cmd_batch_convert()

            assert excinfo.value.code == 1

    @patch("epyr.fair.convert_bruker_to_fair")
    def test_cmd_batch_convert_empty_dir(self, mock_convert, caplog):
        """Test batch convert with empty directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch("sys.argv", ["epyr-batch-convert", temp_dir]):
                with pytest.raises(SystemExit) as excinfo:
                    cli.cmd_batch_convert()

                assert excinfo.value.code == 1

            # Check that error message was logged
            assert any(
                "No .dsc or .spc files found" in record.message
                for record in caplog.records
            )


class TestBaselineCommand:
    """Test baseline correction CLI command."""

    def test_cmd_baseline_missing_file(self, capsys):
        """Test baseline command with missing file."""
        with patch("sys.argv", ["epyr-baseline", "nonexistent.dsc"]):
            with pytest.raises(SystemExit) as excinfo:
                cli.cmd_baseline()

            assert excinfo.value.code == 1

    @patch("epyr.eprload")
    @patch("epyr.baseline.baseline_polynomial")
    def test_cmd_baseline_success(self, mock_baseline, mock_eprload, capsys):
        """Test successful baseline correction."""
        import numpy as np

        # Mock data loading
        x_data = np.linspace(3400, 3500, 100)
        y_data = np.random.randn(100)
        mock_eprload.return_value = (x_data, y_data, {}, "test.dsc")

        # Mock baseline correction
        corrected_y = y_data - 0.1 * x_data
        baseline = 0.1 * x_data
        mock_baseline.return_value = (corrected_y, baseline)

        # Create temporary test file
        with tempfile.NamedTemporaryFile(suffix=".dsc", delete=False) as f:
            test_file = f.name
            f.write(b"test data")

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                output_file = Path(temp_dir) / "test_baseline.csv"

                with patch(
                    "sys.argv",
                    [
                        "epyr-baseline",
                        test_file,
                        "--output",
                        str(output_file),
                        "--method",
                        "polynomial",
                        "--order",
                        "1",
                    ],
                ):
                    cli.cmd_baseline()

                # Check baseline function was called
                mock_baseline.assert_called_once()

                # Check output file was created
                assert output_file.exists()
        finally:
            os.unlink(test_file)
