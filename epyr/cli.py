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

from .config import config
from .eprload import eprload
from .logging_config import get_logger

logger = get_logger(__name__)


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
                       choices=['polynomial', 'exponential', 'stretched_exponential'],
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
            from .baseline import baseline_polynomial
            exclude_regions = args.exclude if args.exclude else None
            y_corrected, baseline = baseline_polynomial(
                y, x_data=x, poly_order=args.order, exclude_regions=exclude_regions
            )
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
    elif args.command == 'validate':
        cmd_validate()
    else:
        parser.print_help()


if __name__ == '__main__':
    main()