"""
Unified Lineshape class for EPR spectroscopy

Provides a single interface for all lineshape types with consistent API.
"""

from typing import Any, Dict, Optional, Tuple, Union

import numpy as np

from .gaussian import gaussian
from .lorentzian import lorentzian
from .lshape import lshape, pseudo_voigt
from .voigtian import voigtian


class Lineshape:
    """
    Unified lineshape class for EPR spectroscopy.

    This class provides a consistent interface for generating all types of
    EPR lineshapes including Gaussian, Lorentzian, Voigt, and pseudo-Voigt
    profiles with support for derivatives, phase rotation, and convolution.

    Parameters:
    -----------
    shape_type : str, default='gaussian'
        Type of lineshape ('gaussian', 'lorentzian', 'voigt', 'pseudo_voigt')
    width : float or tuple
        Full width at half maximum
        - For single shapes: float
        - For Voigt: (gaussian_width, lorentzian_width)
        - For pseudo-Voigt: float (both components use same width)
    alpha : float, default=1.0
        Shape parameter for pseudo-Voigt (1=Gaussian, 0=Lorentzian)
    derivative : int, default=0
        Derivative order (0, 1, 2, or -1 for integral)
    phase : float, default=0.0
        Phase rotation in radians (0=absorption, π/2=dispersion)
    normalize : bool, default=True
        Whether to maintain area normalization

    Examples:
    ---------
    >>> # Create Gaussian lineshape
    >>> gauss = Lineshape('gaussian', width=5.0)
    >>> x = np.linspace(-10, 10, 100)
    >>> y = gauss(x, center=0)
    >>>
    >>> # Lorentzian with first derivative
    >>> lorentz_deriv = Lineshape('lorentzian', width=4.0, derivative=1)
    >>> y_deriv = lorentz_deriv(x, center=2)
    >>>
    >>> # Voigt profile with different widths
    >>> voigt = Lineshape('voigt', width=(3.0, 2.0))
    >>> y_voigt = voigt(x, center=0)
    >>>
    >>> # Pseudo-Voigt 50/50 mix
    >>> pv = Lineshape('pseudo_voigt', width=5.0, alpha=0.5)
    >>> y_pv = pv(x, center=0)
    """

    # Supported lineshape types
    SUPPORTED_SHAPES = {
        "gaussian": gaussian,
        "lorentzian": lorentzian,
        "voigt": voigtian,
        "pseudo_voigt": pseudo_voigt,
        "general": lshape,  # Most flexible option
    }

    def __init__(
        self,
        shape_type: str = "gaussian",
        width: Union[float, Tuple[float, float]] = 1.0,
        alpha: float = 1.0,
        derivative: int = 0,
        phase: float = 0.0,
        normalize: bool = True,
    ):

        # Validate inputs
        if shape_type not in self.SUPPORTED_SHAPES:
            raise ValueError(
                f"shape_type must be one of {list(self.SUPPORTED_SHAPES.keys())}"
            )

        self.shape_type = shape_type
        self.width = width
        self.alpha = alpha
        self.derivative = derivative
        self.phase = phase
        self.normalize = normalize

        # Get the underlying function
        self._func = self.SUPPORTED_SHAPES[shape_type]

        # Store parameters for repr
        self._params = {
            "width": width,
            "alpha": alpha,
            "derivative": derivative,
            "phase": phase,
            "normalize": normalize,
        }

    def __call__(self, x: np.ndarray, center: float, **override_params) -> np.ndarray:
        """
        Generate lineshape at specified points.

        Parameters:
        -----------
        x : array
            Abscissa points
        center : float
            Peak center position
        **override_params : dict
            Parameters to override for this call only

        Returns:
        --------
        array
            Lineshape values
        """
        # Merge default and override parameters
        params = {**self._params, **override_params}

        # Call appropriate function based on shape type
        if self.shape_type == "gaussian":
            return self._func(
                x,
                center,
                params["width"],
                derivative=params["derivative"],
                phase=params["phase"],
            )

        elif self.shape_type == "lorentzian":
            return self._func(
                x,
                center,
                params["width"],
                derivative=params["derivative"],
                phase=params["phase"],
            )

        elif self.shape_type == "voigt":
            return self._func(
                x,
                center,
                params["width"],
                derivative=params["derivative"],
                phase=params["phase"],
            )

        elif self.shape_type in ["pseudo_voigt", "general"]:
            # pseudo_voigt function doesn't support derivative parameter
            if "alpha" in params:
                return self._func(x, center, params["width"], eta=params["alpha"])
            else:
                return self._func(x, center, params["width"])

    def absorption(self, x: np.ndarray, center: float) -> np.ndarray:
        """Generate pure absorption lineshape"""
        return self(x, center, phase=0.0)

    def dispersion(self, x: np.ndarray, center: float) -> np.ndarray:
        """Generate pure dispersion lineshape"""
        return self(x, center, phase=np.pi / 2)

    def derivative(self, x: np.ndarray, center: float, order: int = 1) -> np.ndarray:
        """Generate derivative lineshape"""
        return self(x, center, derivative=order)

    def both_components(
        self, x: np.ndarray, center: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate both absorption and dispersion components.

        Returns:
        --------
        tuple
            (absorption, dispersion) arrays
        """
        if hasattr(self._func, "return_both"):
            # Functions that support return_both parameter
            if self.shape_type in ["gaussian", "lorentzian", "voigt"]:
                return self._func(
                    x, center, self.width, derivative=self.derivative, return_both=True
                )

        # Fallback: compute separately
        abs_part = self.absorption(x, center)
        disp_part = self.dispersion(x, center)
        return abs_part, disp_part

    def set_width(self, width: Union[float, Tuple[float, float]]) -> "Lineshape":
        """Create new Lineshape with different width"""
        return Lineshape(
            self.shape_type,
            width,
            self.alpha,
            self.derivative,
            self.phase,
            self.normalize,
        )

    def set_alpha(self, alpha: float) -> "Lineshape":
        """Create new Lineshape with different alpha (for pseudo-Voigt)"""
        return Lineshape(
            self.shape_type,
            self.width,
            alpha,
            self.derivative,
            self.phase,
            self.normalize,
        )

    def set_derivative(self, derivative: int) -> "Lineshape":
        """Create new Lineshape with different derivative order"""
        return Lineshape(
            self.shape_type,
            self.width,
            self.alpha,
            derivative,
            self.phase,
            self.normalize,
        )

    def set_phase(self, phase: float) -> "Lineshape":
        """Create new Lineshape with different phase"""
        return Lineshape(
            self.shape_type,
            self.width,
            self.alpha,
            self.derivative,
            phase,
            self.normalize,
        )

    def info(self) -> Dict[str, Any]:
        """Get lineshape information"""
        return {
            "shape_type": self.shape_type,
            "width": self.width,
            "alpha": self.alpha,
            "derivative": self.derivative,
            "phase": self.phase,
            "phase_degrees": np.degrees(self.phase),
            "normalize": self.normalize,
            "is_absorption": np.mod(self.phase, 2 * np.pi) == 0,
            "is_dispersion": np.mod(self.phase - np.pi / 2, 2 * np.pi) == 0,
        }

    def __repr__(self) -> str:
        """String representation"""
        return f"Lineshape(type='{self.shape_type}', width={self.width}, α={self.alpha:.2f})"

    def __str__(self) -> str:
        """Human-readable string"""
        phase_desc = "absorption" if self.phase == 0 else f"phase={self.phase:.3f}"
        deriv_desc = "" if self.derivative == 0 else f", d^{self.derivative}"
        return f"{self.shape_type.title()} lineshape (w={self.width}, {phase_desc}{deriv_desc})"


# Factory functions for common lineshapes
def create_gaussian(width: float = 1.0, **kwargs) -> Lineshape:
    """Create Gaussian lineshape"""
    return Lineshape("gaussian", width=width, **kwargs)


def create_lorentzian(width: float = 1.0, **kwargs) -> Lineshape:
    """Create Lorentzian lineshape"""
    return Lineshape("lorentzian", width=width, **kwargs)


def create_voigt(gaussian_width: float, lorentzian_width: float, **kwargs) -> Lineshape:
    """Create Voigt lineshape"""
    return Lineshape("voigt", width=(gaussian_width, lorentzian_width), **kwargs)


def create_pseudo_voigt(width: float = 1.0, alpha: float = 0.5, **kwargs) -> Lineshape:
    """Create pseudo-Voigt lineshape"""
    return Lineshape("pseudo_voigt", width=width, alpha=alpha, **kwargs)
