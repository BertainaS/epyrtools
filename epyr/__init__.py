"""EPyR Tools - Electron Paramagnetic Resonance Tools in Python."""

# Submodule re-exports
from . import baseline
from . import eprplot as plot  # Backward-compatible alias
from . import lineshapes, signalprocessing

# Baseline correction API
from .baseline import (
    MODEL_INFO,
    RegionSelector,
    auto_baseline_with_recommendations,
    baseline_auto_1d,
    baseline_bi_exponential_1d,
    baseline_polynomial_1d,
    baseline_polynomial_2d,
    baseline_stretched_exponential_1d,
    bi_exponential_1d,
    close_selector_window,
    compare_models_detailed,
    configure,
    create_region_mask_1d,
    create_region_mask_2d,
    exponential_1d,
    get_baseline_regions_1d,
    get_baseline_regions_2d,
    get_configuration,
    get_model_function,
    get_model_recommendations,
    interactive_select_regions_1d,
    interactive_select_regions_2d,
    is_interactive_available,
    jupyter_help,
    list_available_models,
    polynomial_1d,
    polynomial_2d,
    setup_inline_backend,
    setup_interactive_backend,
    setup_notebook_backend,
    setup_widget_backend,
    stretched_exponential_1d,
    validate_regions_1d,
    validate_regions_2d,
)
from .config import config
from .eprload import eprload

# Plotting API
from .eprplot import plot_1d, plot_2d_map, plot_2d_slicer, plot_2d_waterfall

# FAIR conversion API
from .fair import (
    BRUKER_PARAM_MAP,
    ValidationResult,
    append_fair_metadata,
    batch_convert_directory,
    convert_bruker_to_fair,
    create_validation_report,
    extract_axis_info,
    process_parameters,
    save_fair,
    save_to_csv,
    save_to_csv_json,
    save_to_hdf5,
    save_to_jpg,
    save_to_json,
    validate_fair_dataset,
)
from .isotope_gui import run_gui as isotopes
from .lineshapes import Lineshape, gaussian, lorentzian, pseudo_voigt, voigtian
from .logging_config import get_logger, setup_logging
from .performance import DataCache, OptimizedLoader, get_performance_info

# Physics / units API
from .physics import (
    AVOGADRO,
    AVOGADRO_CGS,
    BMAGN,
    BMAGN_CGS,
    BOLTZM,
    BOLTZM_CGS,
    CLIGHT,
    CLIGHT_CGS,
    ECHARGE,
    ECHARGE_CGS,
    EVOLT,
    EVOLT_CGS,
    GFREE,
    GFREE_CGS,
    HBAR,
    HBAR_CGS,
    NMAGN,
    NMAGN_CGS,
    PLANCK,
    PLANCK_CGS,
    avogadro,
    bmagn,
    boltzm,
    clight,
    cm_inv_to_mhz,
    constants,
    constants_summary,
    demo_conversions,
    echarge,
    energy_conversion_table,
    evolt,
    frequency_field_conversion_table,
    frequency_to_magnetic_field,
    gamma_hz,
    gfree,
    hbar,
    list_conversions,
    magnetic_field_to_frequency,
    mhz_to_cm_inv,
    mhz_to_mt,
    mt_to_mhz,
    nmagn,
    planck,
    thermal_energy,
    unitconvert,
    wavelength_to_frequency,
)
from .plugins import plugin_manager
from .sub.utils import BrukerListFiles

__version__ = "0.3.9"

logger = get_logger(__name__)
logger.debug("Package 'epyr' v%s initialized.", __version__)
