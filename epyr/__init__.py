"""EPyR Tools - Electron Paramagnetic Resonance Tools in Python."""

# Import configuration and logging first
# Import baseline module for compatibility and convenience
from . import baseline

# Import backend control functions for convenience
# Import baseline correction functions from new modular package
from .baseline import *
from .baseline import setup_inline_backend, setup_notebook_backend, setup_widget_backend
from .config import config
from .logging_config import get_logger, setup_logging

from . import eprplot as plot  # Alias for backward compatibility
from . import lineshapes, signalprocessing
from .eprload import *
from .eprplot import *
from .fair import *
from .isotope_gui import run_gui as isotopes
from .lineshapes import Lineshape, gaussian, lorentzian, pseudo_voigt, voigtian
from .performance import DataCache, OptimizedLoader, get_performance_info
from .physics import *
from .physics import constants
from .plugins import plugin_manager
from .sub.utils import BrukerListFiles

__version__ = "0.3.5"

# Set up logging
logger = get_logger(__name__)

# Display version on import
print(f"EPyR Tools v{__version__}")
logger.debug(f"Package 'epyr' v{__version__} initialized.")
