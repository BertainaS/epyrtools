"""EPyR Tools - Electron Paramagnetic Resonance Tools in Python."""

# Import configuration and logging first
# Import baseline module for compatibility and convenience
from . import baseline
from . import eprplot as plot  # Alias for backward compatibility
from . import lineshapes, signalprocessing

# Import backend control functions for convenience
# Import baseline correction functions from new modular package
from .baseline import *
from .baseline import setup_inline_backend, setup_notebook_backend, setup_widget_backend
from .config import config
from .eprload import eprload
from .eprplot import *
from .fair import *
from .isotope_gui import run_gui as isotopes
from .lineshapes import Lineshape, gaussian, lorentzian, pseudo_voigt, voigtian
from .logging_config import get_logger, setup_logging
from .performance import DataCache, OptimizedLoader, get_performance_info
from .physics import *
from .physics import constants
from .plugins import plugin_manager
from .sub.utils import BrukerListFiles

__version__ = "0.3.8"

# Set up logging
logger = get_logger(__name__)

logger.debug("Package 'epyr' v%s initialized.", __version__)
