"""EPyR Tools - Electron Paramagnetic Resonance Tools in Python."""

# Import configuration and logging first
from .config import config
from .logging_config import setup_logging, get_logger

# Import specific, useful components from the modules
from . import baseline
from . import lineshapes
from .constants import *
from .eprload import *
from .fair import *
from .isotope_gui import run_gui as isotopes
from .lineshapes import Lineshape, gaussian, lorentzian, voigtian, pseudo_voigt
from .performance import OptimizedLoader, DataCache, get_performance_info
from .plugins import plugin_manager
from .plot import *
from .sub.utils import BrukerListFiles

__version__ = "0.1.6"

# Set up logging
logger = get_logger(__name__)
logger.debug(f"Package 'epyr' v{__version__} initialized.")
