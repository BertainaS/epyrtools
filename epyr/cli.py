"""
Command Line Interface for EPyR Tools
=====================================

Provides command-line tools for common EPyR workflows:
- Data conversion (Bruker -> FAIR formats)
- Baseline correction
- Batch processing
- Configuration management

Usage:
    epyr-convert input.dsc --output-dir ./results
    epyr-baseline spectrum.dsc --method polynomial --order 2
    epyr-batch-convert ./data/ --format csv,json
    epyr-config --set plotting.dpi 300
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

from .config import config
from .eprload import eprload
from .logging_config import get_logger

logger = get_logger(__name__)


class InteractiveMeasurementTool:
    """Interactive tool for measuring distances between two points on a plot."""
    
    def __init__(self, ax, x_data, y_data):
        """Initialize the measurement tool.
        
        Args:
            ax: Matplotlib axes object
            x_data: X-axis data array
            y_data: Y-axis data array
        """
        self.ax = ax
        self.x_data = x_data
        self.y_data = y_data
        self.points = []
        self.lines = []
        self.annotations = []
        self.cid = None
        
    def enable(self):
        """Enable the measurement tool."""
        self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_click)
        print("\n📏 Measurement tool enabled!")
        print("Instructions:")
        print("  • Click two points on the plot to measure distance")
        print("  • Right-click to clear measurements")
        print("  • Press 'q' or close window to exit")
        
    def disable(self):
        """Disable the measurement tool."""
        if self.cid:
            self.ax.figure.canvas.mpl_disconnect(self.cid)
            self.cid = None
            
    def on_click(self, event):
        """Handle mouse click events."""
        if event.inaxes != self.ax:
            return
            
        if event.button == 3:  # Right click - clear measurements
            self.clear_measurements()
            return
            
        if event.button != 1:  # Only handle left clicks
            return
            
        # Add point
        x, y = event.xdata, event.ydata
        if x is None or y is None:
            return
            
        self.points.append((x, y))
        
        # Plot the point
        point_plot = self.ax.plot(x, y, 'ro', markersize=8, markeredgecolor='white', markeredgewidth=1)[0]
        self.lines.append(point_plot)
        
        print(f"Point {len(self.points)}: x={x:.4f}, y={y:.4e}")
        
        # If we have two points, calculate and display distance
        if len(self.points) == 2:
            self.calculate_distance()
            self.points = []  # Reset for next measurement
            
        self.ax.figure.canvas.draw()
        
    def calculate_distance(self):
        """Calculate and display distance between two points."""
        if len(self.points) != 2:
            return
            
        p1, p2 = self.points
        x1, y1 = p1
        x2, y2 = p2
        
        # Calculate deltas
        delta_x = x2 - x1
        delta_y = y2 - y1
        distance = np.sqrt(delta_x**2 + delta_y**2)
        
        # Draw line between points
        line = self.ax.plot([x1, x2], [y1, y2], 'r--', linewidth=2, alpha=0.7)[0]
        self.lines.append(line)
        
        # Add measurement annotation
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        
        measurement_text = f'Δx = {delta_x:.4f}\nΔy = {delta_y:.4e}\n|Δ| = {distance:.4e}'
        annotation = self.ax.annotate(
            measurement_text,
            xy=(mid_x, mid_y),
            xytext=(10, 10),
            textcoords='offset points',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8),
            fontsize=9,
            ha='left'
        )
        self.annotations.append(annotation)
        
        # Print results to console
        print(f"\n📐 Measurement Results:")
        print(f"  Point 1: ({x1:.4f}, {y1:.4e})")
        print(f"  Point 2: ({x2:.4f}, {y2:.4e})")
        print(f"  Δx = {delta_x:.4f}")
        print(f"  Δy = {delta_y:.4e}")
        print(f"  Distance = {distance:.4e}")
        print(f"\nClick two more points for another measurement, or right-click to clear.\n")
        
    def clear_measurements(self):
        """Clear all measurements from the plot."""
        # Remove all plotted elements
        for item in self.lines + self.annotations:
            if item in self.ax.lines:
                item.remove()
            elif item in self.ax.texts:
                item.remove()
                
        self.lines.clear()
        self.annotations.clear()
        self.points.clear()
        
        self.ax.figure.canvas.draw()
        print("🧹 Measurements cleared. Click two points for a new measurement.")


def create_interactive_plot_with_measurements(x, y, params, file_path, enable_measurements=False):
    """Create an interactive plot with optional measurement tools.
    
    Args:
        x: X-axis data
        y: Y-axis data  
        params: Parameter dictionary
        file_path: Path to the loaded file
        enable_measurements: Whether to enable measurement tool
    """
    import matplotlib.pyplot as plt
    from pathlib import Path
    
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_title = Path(file_path).name if file_path else "EPR Data"
    ax.set_title(plot_title, fontsize=12)
    
    # Plot the data
    if y.ndim == 1:
        # 1D data
        absc = x if x is not None and hasattr(x, '__len__') else np.arange(len(y))
        
        if np.isrealobj(y):
            ax.plot(absc, y, 'b-', linewidth=1.5, label='data')
        else:
            ax.plot(absc, np.real(y), 'b-', linewidth=1.5, label='real')
            ax.plot(absc, np.imag(y), 'r--', linewidth=1.5, label='imag')
            ax.legend()
            
        # Set labels
        x_label = params.get("XAXIS_NAME", "Field") if params else "Field"
        x_unit = params.get("XAXIS_UNIT", "G") if params else "G"
        if x_unit:
            x_label += f" ({x_unit})"
            
        ax.set_xlabel(x_label)
        ax.set_ylabel("Intensity (a.u.)")
        ax.grid(True, linestyle=":", alpha=0.6)
        
        # Enable measurement tool if requested
        measurement_tool = None
        if enable_measurements:
            measurement_tool = InteractiveMeasurementTool(ax, absc, y)
            measurement_tool.enable()
            
    else:
        # 2D data - basic implementation
        ax.imshow(np.real(y), aspect='auto', cmap='viridis')
        ax.set_title(f"{plot_title} (2D data)")
        print("📊 2D data plotted. Measurement tool works best with 1D data.")
        
    plt.tight_layout()
    
    # Add keyboard shortcuts
    def on_key(event):
        if event.key == 'q':
            plt.close('all')
        elif event.key == 'c' and enable_measurements and measurement_tool:
            measurement_tool.clear_measurements()
            
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    if enable_measurements:
        print("\n⌨️  Keyboard shortcuts:")
        print("  • 'c' - Clear measurements")
        print("  • 'q' - Quit")
    
    return fig, ax


def cmd_convert():
    """Convert Bruker files to FAIR formats."""
    parser = argparse.ArgumentParser(
        prog='epyr-convert',
        description='Convert Bruker EPR files to FAIR formats (CSV, JSON, HDF5)'
    )
    parser.add_argument('input', help='Input Bruker file (.dta, .dsc, .spc, .par)')
    parser.add_argument('-o', '--output-dir', default='.', 
                       help='Output directory (default: current directory)')
    parser.add_argument('-f', '--formats', default='csv,json',
                       help='Output formats: csv,json,hdf5 (default: csv,json)')
    parser.add_argument('--no-metadata', action='store_true',
                       help='Skip metadata export')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    formats = [f.strip().lower() for f in args.formats.split(',')]
    
    try:
        from .fair import convert_bruker_to_fair
        
        logger.info(f"Converting {input_path} to formats: {', '.join(formats)}")
        
        # Convert to specified formats
        success = convert_bruker_to_fair(
            str(input_path),
            output_dir=str(output_dir),
            formats=formats,
            include_metadata=not args.no_metadata
        )
        
        if success:
            logger.info(f"Conversion completed successfully. Output in: {output_dir}")
        else:
            logger.error("Conversion failed")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Conversion error: {e}")
        if args.verbose:
            logger.debug("Full traceback:", exc_info=True)
        sys.exit(1)


def cmd_baseline():
    """Apply baseline correction to EPR data."""
    parser = argparse.ArgumentParser(
        prog='epyr-baseline',
        description='Apply baseline correction to EPR data'
    )
    parser.add_argument('input', help='Input EPR file')
    parser.add_argument('-o', '--output', help='Output file (default: input_baseline.csv)')
    parser.add_argument('-m', '--method', default='polynomial', 
                       choices=['polynomial', 'stretched_exponential', 'bi_exponential', 'auto'],
                       help='Baseline correction method')
    parser.add_argument('--order', type=int, default=1,
                       help='Polynomial order (for polynomial method)')
    parser.add_argument('--exclude', action='append', nargs=2, type=float,
                       metavar=('START', 'END'),
                       help='Exclude region from fit (can be used multiple times)')
    parser.add_argument('--plot', action='store_true',
                       help='Generate comparison plot')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(f"{input_path.stem}_baseline.csv")
    
    try:
        # Load data
        logger.info(f"Loading data from {input_path}")
        x, y, params, _ = eprload(str(input_path), plot_if_possible=False)
        
        if x is None or y is None:
            logger.error("Failed to load data")
            sys.exit(1)
        
        # Apply baseline correction
        logger.info(f"Applying {args.method} baseline correction")
        
        if args.method == 'polynomial':
            from .baseline import baseline_polynomial_1d
            exclude_regions = args.exclude if args.exclude else None
            # Convert exclude_regions format for new API
            manual_regions = exclude_regions
            region_mode = 'exclude' if manual_regions else None
            
            y_corrected, baseline = baseline_polynomial_1d(
                x, y, params, 
                order=args.order, 
                manual_regions=manual_regions,
                region_mode=region_mode
            )
        elif args.method == 'stretched_exponential':
            from .baseline import baseline_stretched_exponential_1d
            y_corrected, baseline = baseline_stretched_exponential_1d(x, y, params)
        elif args.method == 'bi_exponential':
            from .baseline import baseline_bi_exponential_1d
            y_corrected, baseline = baseline_bi_exponential_1d(x, y, params)
        elif args.method == 'auto':
            from .baseline import baseline_auto_1d
            y_corrected, baseline, info = baseline_auto_1d(x, y, params, verbose=True)
            logger.info(f"Automatic selection chose: {info['best_model']}")
        else:
            logger.error(f"Method {args.method} not yet implemented in CLI")
            sys.exit(1)
        
        # Save results
        import pandas as pd
        df = pd.DataFrame({
            'field': x if hasattr(x, '__len__') else range(len(y)),
            'original': y,
            'baseline': baseline,
            'corrected': y_corrected
        })
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        # Generate plot if requested
        if args.plot:
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            field = x if hasattr(x, '__len__') else range(len(y))
            
            ax1.plot(field, y, 'b-', label='Original', alpha=0.7)
            ax1.plot(field, baseline, 'r--', label='Baseline')
            ax1.plot(field, y_corrected, 'g-', label='Corrected')
            ax1.legend()
            ax1.set_title('Baseline Correction')
            ax1.grid(True, alpha=0.3)
            
            ax2.plot(field, y_corrected, 'g-', linewidth=2)
            ax2.set_title('Corrected Spectrum')
            ax2.set_xlabel('Field' if hasattr(x, '__len__') else 'Index')
            ax2.set_ylabel('Intensity')
            ax2.grid(True, alpha=0.3)
            
            plot_path = output_path.with_suffix('.png')
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300)
            logger.info(f"Plot saved to {plot_path}")
            
    except Exception as e:
        logger.error(f"Baseline correction error: {e}")
        if args.verbose:
            logger.debug("Full traceback:", exc_info=True)
        sys.exit(1)


def cmd_batch_convert():
    """Batch convert multiple files."""
    parser = argparse.ArgumentParser(
        prog='epyr-batch-convert',
        description='Batch convert multiple Bruker EPR files'
    )
    parser.add_argument('input_dir', help='Input directory containing Bruker files')
    parser.add_argument('-o', '--output-dir', help='Output directory (default: input_dir/converted)')
    parser.add_argument('-f', '--formats', default='csv,json',
                       help='Output formats: csv,json,hdf5')
    parser.add_argument('--pattern', default='*.dsc',
                       help='File pattern to match (default: *.dsc)')
    parser.add_argument('-j', '--jobs', type=int, default=1,
                       help='Number of parallel jobs (default: 1)')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)
    
    output_dir = Path(args.output_dir) if args.output_dir else input_dir / 'converted'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find files to convert
    files = list(input_dir.glob(args.pattern))
    if not files:
        logger.error(f"No files found matching pattern: {args.pattern}")
        sys.exit(1)
    
    logger.info(f"Found {len(files)} files to convert")
    
    formats = [f.strip().lower() for f in args.formats.split(',')]
    
    # Convert files
    success_count = 0
    for i, file_path in enumerate(files, 1):
        try:
            logger.info(f"[{i}/{len(files)}] Converting {file_path.name}")
            
            from .fair import convert_bruker_to_fair
            success = convert_bruker_to_fair(
                str(file_path),
                output_dir=str(output_dir),
                formats=formats
            )
            
            if success:
                success_count += 1
            else:
                logger.warning(f"Failed to convert {file_path.name}")
                
        except Exception as e:
            logger.error(f"Error converting {file_path.name}: {e}")
            if args.verbose:
                logger.debug("Full traceback:", exc_info=True)
    
    logger.info(f"Batch conversion completed: {success_count}/{len(files)} files converted")


def cmd_config():
    """Configuration management."""
    parser = argparse.ArgumentParser(
        prog='epyr-config',
        description='Manage EPyR Tools configuration'
    )
    
    subparsers = parser.add_subparsers(dest='action', help='Configuration actions')
    
    # Show config
    show_parser = subparsers.add_parser('show', help='Show current configuration')
    show_parser.add_argument('section', nargs='?', help='Configuration section to show')
    
    # Set config
    set_parser = subparsers.add_parser('set', help='Set configuration value')
    set_parser.add_argument('key', help='Configuration key (e.g., plotting.dpi)')
    set_parser.add_argument('value', help='Configuration value')
    
    # Reset config
    reset_parser = subparsers.add_parser('reset', help='Reset configuration')
    reset_parser.add_argument('section', nargs='?', help='Section to reset (or all)')
    
    # Export/Import
    export_parser = subparsers.add_parser('export', help='Export configuration')
    export_parser.add_argument('file', help='Output file')
    
    import_parser = subparsers.add_parser('import', help='Import configuration')
    import_parser.add_argument('file', help='Input file')
    
    args = parser.parse_args()
    
    if not args.action:
        parser.print_help()
        return
    
    try:
        if args.action == 'show':
            if args.section:
                section_config = config.get_section(args.section)
                if section_config:
                    import json
                    print(json.dumps(section_config, indent=2))
                else:
                    logger.error(f"Section '{args.section}' not found")
            else:
                import json
                print(json.dumps(config._config, indent=2))
        
        elif args.action == 'set':
            # Try to parse value as JSON first
            try:
                import json
                value = json.loads(args.value)
            except json.JSONDecodeError:
                value = args.value
            
            config.set(args.key, value)
            config.save()
            logger.info(f"Set {args.key} = {value}")
        
        elif args.action == 'reset':
            if args.section and args.section != 'all':
                config.reset_section(args.section)
                logger.info(f"Reset section: {args.section}")
            else:
                config.reset_all()
                logger.info("Reset all configuration to defaults")
            config.save()
        
        elif args.action == 'export':
            config.export_config(args.file)
            logger.info(f"Configuration exported to {args.file}")
        
        elif args.action == 'import':
            config.import_config(args.file)
            config.save()
            logger.info(f"Configuration imported from {args.file}")
        
    except Exception as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)


def cmd_info():
    """Show system and configuration information."""
    parser = argparse.ArgumentParser(
        prog='epyr-info',
        description='Display EPyR Tools system and configuration information'
    )
    parser.add_argument('--config', action='store_true',
                       help='Show configuration details')
    parser.add_argument('--performance', action='store_true', 
                       help='Show performance information')
    parser.add_argument('--plugins', action='store_true',
                       help='Show loaded plugins')
    parser.add_argument('--all', action='store_true',
                       help='Show all information')
    
    args = parser.parse_args()
    
    import json
    from . import __version__
    
    # Show version info
    print(f"EPyR Tools Version: {__version__}")
    print(f"Configuration file: {config.get_config_file_path()}")
    print()
    
    if args.config or args.all:
        print("=== Configuration ===")
        print(json.dumps(config._config, indent=2))
        print()
    
    if args.performance or args.all:
        print("=== Performance Information ===")
        from .performance import get_performance_info
        perf_info = get_performance_info()
        print(json.dumps(perf_info, indent=2))
        print()
    
    if args.plugins or args.all:
        print("=== Loaded Plugins ===")
        from .plugins import plugin_manager
        plugins_info = plugin_manager.list_plugins()
        print(json.dumps(plugins_info, indent=2))
        print()


def cmd_isotopes():
    """Launch the isotope database GUI."""
    parser = argparse.ArgumentParser(
        prog='epyr-isotopes',
        description='Launch the interactive isotope database GUI'
    )
    
    args = parser.parse_args()
    
    try:
        logger.info("Launching isotope database GUI...")
        from .isotope_gui import run_gui
        run_gui()
    except Exception as e:
        logger.error(f"Failed to launch isotope GUI: {e}")
        sys.exit(1)


def _plot_main(args_list=None):
    """Main plotting function that can accept custom args."""
    parser = argparse.ArgumentParser(
        prog='epyr-plot',
        description='Load and plot EPR data files with interactive visualization'
    )
    parser.add_argument('file', nargs='?', help='EPR file to load (.dta, .dsc, .spc, .par). If not provided, opens file dialog.')
    parser.add_argument('-s', '--scaling', default='', 
                       help='Scaling string (n=scans, P=power, G=gain, T=temp, c=time)')
    parser.add_argument('--no-plot', action='store_true',
                       help='Load data without plotting')
    parser.add_argument('--interactive', action='store_true',
                       help='Enable interactive matplotlib backend')
    parser.add_argument('--save', action='store_true',
                       help='Save plot as PNG file')
    parser.add_argument('--measure', action='store_true',
                       help='Enable interactive measurement tool (click two points to measure distance)')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args(args_list)
    
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    # Set up interactive backend if requested
    if args.interactive:
        import matplotlib
        import platform
        
        if platform.system() == "Darwin":  # macOS
            try:
                matplotlib.use('TkAgg')
                logger.info("Using TkAgg backend for interactive plotting on macOS")
            except ImportError:
                logger.warning("TkAgg not available, using default backend")
        else:
            try:
                matplotlib.use('Qt5Agg')
                logger.info("Using Qt5Agg backend for interactive plotting")
            except ImportError:
                try:
                    matplotlib.use('TkAgg')
                    logger.info("Using TkAgg backend for interactive plotting")
                except ImportError:
                    logger.warning("No interactive backend available, using default")
    
    try:
        # Load the data
        logger.info("Loading EPR data...")
        
        # For measurement mode, disable default plotting
        plot_with_eprload = not args.no_plot and not (args.interactive and args.measure)
        
        x, y, params, file_path = eprload(
            args.file,
            scaling=args.scaling,
            plot_if_possible=plot_with_eprload,
            save_if_possible=args.save and not args.measure
        )
        
        if x is None or y is None:
            logger.error("Failed to load data or loading was cancelled")
            sys.exit(1)
        
        logger.info(f"Successfully loaded: {file_path}")
        logger.info(f"Data shape: {y.shape}")
        
        if hasattr(x, 'shape'):
            logger.info(f"X-axis shape: {x.shape}")
        elif isinstance(x, (list, tuple)):
            logger.info(f"X-axis shapes: {[ax.shape for ax in x]}")
        
        logger.info(f"Parameters loaded: {len(params) if params else 0}")
        
        # Show key parameters
        if params:
            key_params = ['MWFQ', 'MWPW', 'RCAG', 'AVGS', 'SPTP']
            found_params = {k: params.get(k) for k in key_params if k in params}
            if found_params:
                logger.info("Key parameters:")
                for k, v in found_params.items():
                    logger.info(f"  {k}: {v}")
        
        # Handle interactive plotting with optional measurements
        if args.interactive and not args.no_plot:
            if args.measure:
                # Use custom interactive plot with measurement tools
                logger.info("Creating interactive plot with measurement tools...")
                fig, ax = create_interactive_plot_with_measurements(
                    x, y, params, file_path, enable_measurements=True
                )
                
                if args.save:
                    from pathlib import Path
                    save_path = Path(file_path).with_suffix('.png') if file_path else Path('epr_plot.png')
                    fig.savefig(save_path, dpi=300)
                    logger.info(f"Plot saved to {save_path}")
                
                import matplotlib.pyplot as plt
                plt.show(block=True)
                logger.info("Interactive measurement plot closed.")
            else:
                # Standard interactive plot
                import matplotlib.pyplot as plt
                plt.show(block=True)
                logger.info("Interactive plot displayed. Close the plot window to exit.")
        
    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        if args.verbose:
            logger.debug("Full traceback:", exc_info=True)
        sys.exit(1)


def cmd_plot():
    """Load and plot EPR data files interactively."""
    _plot_main()


def cmd_plot_with_args(args):
    """Load and plot EPR data files interactively with pre-parsed args."""
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    # Set up interactive backend if requested
    if args.interactive:
        import matplotlib
        import platform
        
        if platform.system() == "Darwin":  # macOS
            try:
                matplotlib.use('TkAgg')
                logger.info("Using TkAgg backend for interactive plotting on macOS")
            except ImportError:
                logger.warning("TkAgg not available, using default backend")
        else:
            try:
                matplotlib.use('Qt5Agg')
                logger.info("Using Qt5Agg backend for interactive plotting")
            except ImportError:
                try:
                    matplotlib.use('TkAgg')
                    logger.info("Using TkAgg backend for interactive plotting")
                except ImportError:
                    logger.warning("No interactive backend available, using default")
    
    try:
        # Load the data
        logger.info("Loading EPR data...")
        
        # For measurement mode, disable default plotting
        plot_with_eprload = not args.no_plot and not (args.interactive and args.measure)
        
        x, y, params, file_path = eprload(
            args.file,
            scaling=args.scaling,
            plot_if_possible=plot_with_eprload,
            save_if_possible=args.save and not args.measure
        )
        
        if x is None or y is None:
            logger.error("Failed to load data or loading was cancelled")
            sys.exit(1)
        
        logger.info(f"Successfully loaded: {file_path}")
        logger.info(f"Data shape: {y.shape}")
        
        if hasattr(x, 'shape'):
            logger.info(f"X-axis shape: {x.shape}")
        elif isinstance(x, (list, tuple)):
            logger.info(f"X-axis shapes: {[ax.shape for ax in x]}")
        
        logger.info(f"Parameters loaded: {len(params) if params else 0}")
        
        # Show key parameters
        if params:
            key_params = ['MWFQ', 'MWPW', 'RCAG', 'AVGS', 'SPTP']
            found_params = {k: params.get(k) for k in key_params if k in params}
            if found_params:
                logger.info("Key parameters:")
                for k, v in found_params.items():
                    logger.info(f"  {k}: {v}")
        
        # Handle interactive plotting with optional measurements
        if args.interactive and not args.no_plot:
            if args.measure:
                # Use custom interactive plot with measurement tools
                logger.info("Creating interactive plot with measurement tools...")
                fig, ax = create_interactive_plot_with_measurements(
                    x, y, params, file_path, enable_measurements=True
                )
                
                if args.save:
                    from pathlib import Path
                    save_path = Path(file_path).with_suffix('.png') if file_path else Path('epr_plot.png')
                    fig.savefig(save_path, dpi=300)
                    logger.info(f"Plot saved to {save_path}")
                
                import matplotlib.pyplot as plt
                plt.show(block=True)
                logger.info("Interactive measurement plot closed.")
            else:
                # Standard interactive plot
                import matplotlib.pyplot as plt
                plt.show(block=True)
                logger.info("Interactive plot displayed. Close the plot window to exit.")
        
    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        if args.verbose:
            logger.debug("Full traceback:", exc_info=True)
        sys.exit(1)


def cmd_validate():
    """Validate EPR data files."""
    parser = argparse.ArgumentParser(
        prog='epyr-validate',
        description='Validate EPR data files for integrity and format compliance'
    )
    parser.add_argument('files', nargs='+', help='Files to validate')
    parser.add_argument('--format', help='Expected file format (auto-detect if not specified)')
    parser.add_argument('--detailed', action='store_true',
                       help='Show detailed validation results')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        from .logging_config import setup_logging
        setup_logging('DEBUG')
    
    total_files = len(args.files)
    valid_files = 0
    
    for file_path in args.files:
        file_path = Path(file_path)
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            continue
        
        try:
            # Try to load the file
            logger.info(f"Validating {file_path}")
            x, y, params, _ = eprload(str(file_path), plot_if_possible=False)
            
            if x is not None and y is not None:
                # Perform FAIR validation if detailed output requested
                if args.detailed:
                    from .fair.validation import validate_fair_dataset
                    
                    data_dict = {
                        'x_data': x,
                        'y_data': y,
                        'metadata': params or {}
                    }
                    
                    fair_result = validate_fair_dataset(data_dict, file_path)
                    
                    if fair_result.is_valid:
                        print(f"✓ {file_path.name} - Valid")
                        valid_files += 1
                    else:
                        print(f"⚠ {file_path.name} - Valid data but FAIR compliance issues")
                        valid_files += 1
                    
                    print(f"  Data points: {len(y)}")
                    print(f"  X-axis range: {np.min(x) if x is not None else 'N/A'} to {np.max(x) if x is not None else 'N/A'}")
                    print(f"  Parameters: {len(params) if params else 0} entries")
                    print(f"  FAIR compliance: {len(fair_result.errors)} errors, {len(fair_result.warnings)} warnings")
                    
                    if fair_result.errors:
                        for error in fair_result.errors[:3]:  # Show first 3 errors
                            print(f"    Error: {error}")
                        if len(fair_result.errors) > 3:
                            print(f"    ... and {len(fair_result.errors) - 3} more errors")
                else:
                    valid_files += 1
                    print(f"✓ {file_path.name} - Valid")
            else:
                logger.warning(f"Failed to extract valid data from {file_path}")
                print(f"✗ {file_path.name} - Invalid data")
                
        except Exception as e:
            logger.error(f"Validation failed for {file_path}: {e}")
            print(f"✗ {file_path.name} - Error: {e}")
    
    print(f"\nValidation Summary: {valid_files}/{total_files} files valid")
    
    if valid_files < total_files:
        sys.exit(1)


def main():
    """Main CLI entry point - shows available commands."""
    parser = argparse.ArgumentParser(
        prog='epyr',
        description='EPyR Tools - Command Line Interface'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add subcommands
    subparsers.add_parser('convert', help='Convert Bruker files to FAIR formats')
    subparsers.add_parser('baseline', help='Apply baseline correction')
    subparsers.add_parser('batch-convert', help='Batch convert multiple files')
    subparsers.add_parser('config', help='Configuration management')
    subparsers.add_parser('info', help='Show system and configuration info')
    subparsers.add_parser('isotopes', help='Launch isotope database GUI')
    
    # Plot subcommand with arguments
    plot_parser = subparsers.add_parser('plot', help='Load and plot EPR data interactively')
    plot_parser.add_argument('file', nargs='?', help='EPR file to load (.dta, .dsc, .spc, .par). If not provided, opens file dialog.')
    plot_parser.add_argument('-s', '--scaling', default='', 
                           help='Scaling string (n=scans, P=power, G=gain, T=temp, c=time)')
    plot_parser.add_argument('--no-plot', action='store_true',
                           help='Load data without plotting')
    plot_parser.add_argument('--interactive', action='store_true',
                           help='Enable interactive matplotlib backend')
    plot_parser.add_argument('--save', action='store_true',
                           help='Save plot as PNG file')
    plot_parser.add_argument('--measure', action='store_true',
                           help='Enable interactive measurement tool (click two points to measure distance)')
    plot_parser.add_argument('-v', '--verbose', action='store_true')
    
    subparsers.add_parser('validate', help='Validate EPR data files')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        print("\nUse 'epyr <command> --help' for more information on a specific command.")
        return
    
    # Dispatch to appropriate command
    if args.command == 'convert':
        cmd_convert()
    elif args.command == 'baseline':
        cmd_baseline()
    elif args.command == 'batch-convert':
        cmd_batch_convert()
    elif args.command == 'config':
        cmd_config()
    elif args.command == 'info':
        cmd_info()
    elif args.command == 'isotopes':
        cmd_isotopes()
    elif args.command == 'plot':
        cmd_plot_with_args(args)
    elif args.command == 'validate':
        cmd_validate()
    else:
        parser.print_help()


if __name__ == '__main__':
    main()