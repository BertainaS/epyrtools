"""
EPyR Tools - Lineshape Functions Module

Comprehensive collection of EPR lineshape functions for spectroscopy analysis.
Includes Gaussian, Lorentzian, Voigt, and pseudo-Voigt profiles with:

- Derivatives (1st and 2nd order)
- Phase rotation (absorption/dispersion modes)
- Spectrum convolution capabilities
- Modern, optimized implementations

Classes and Functions:
    Lineshape: Main lineshape class with all profiles
    gaussian: Pure Gaussian profiles
    lorentzian: Pure Lorentzian profiles
    voigtian: True Voigt profiles (convolution)
    pseudo_voigt: Pseudo-Voigt (linear combination)
    convspec: Spectrum convolution

Author: EPyR Tools Development Team
License: MIT
"""

from .convspec import convspec
from .gaussian import (
    gaussian,
    gaussian_absorption,
    gaussian_derivative,
    gaussian_dispersion,
)
from .lineshape_class import (
    Lineshape,
    create_gaussian,
    create_lorentzian,
    create_pseudo_voigt,
    create_voigt,
)
from .lorentzian import (
    lorentzian,
    lorentzian_absorption,
    lorentzian_derivative,
    lorentzian_dispersion,
)
from .lshape import lshape, pseudo_voigt
from .voigtian import voigtian

# Version info
__version__ = "0.3.1"
__author__ = "EPyR Tools Development Team"

# Import fitting functionality
from .fitting import FitResult, fit_epr_signal, fit_multiple_shapes

# Main functions for easy access
__all__ = [
    # Main class
    "Lineshape",
    # Individual lineshape functions
    "gaussian",
    "lorentzian",
    "voigtian",
    "lshape",
    "pseudo_voigt",
    "convspec",
    # Factory functions
    "create_gaussian",
    "create_lorentzian",
    "create_voigt",
    "create_pseudo_voigt",
    # Convenience functions
    "gaussian_absorption",
    "gaussian_dispersion",
    "gaussian_derivative",
    "lorentzian_absorption",
    "lorentzian_dispersion",
    "lorentzian_derivative",
    # Fitting functionality
    "fit_epr_signal",
    "fit_multiple_shapes",
    "FitResult",
]


# Module-level convenience function
def create_lineshape(shape_type="gaussian", **kwargs):
    """
    Create a lineshape function with specified type.

    Parameters:
    -----------
    shape_type : str
        Type of lineshape ('gaussian', 'lorentzian', 'voigt', 'pseudo_voigt')
    **kwargs :
        Arguments passed to lineshape function

    Returns:
    --------
    function
        Configured lineshape function

    Examples:
    ---------
    >>> # Create a Gaussian lineshape
    >>> gauss_func = create_lineshape('gaussian', width=5.0)
    >>> x = np.linspace(-10, 10, 100)
    >>> y = gauss_func(x, center=0)
    """
    shape_map = {
        "gaussian": gaussian,
        "lorentzian": lorentzian,
        "voigt": voigtian,
        "pseudo_voigt": pseudo_voigt,
        "general": lshape,
    }

    if shape_type not in shape_map:
        raise ValueError(
            f"Unknown shape_type: {shape_type}. Choose from {list(shape_map.keys())}"
        )

    shape_func = shape_map[shape_type]

    # Return a partial function with preset kwargs
    def configured_lineshape(x, center, **extra_kwargs):
        combined_kwargs = {**kwargs, **extra_kwargs}
        return shape_func(x, center, **combined_kwargs)

    return configured_lineshape
