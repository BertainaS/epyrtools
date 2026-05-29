"""Mathematical models for baseline correction.

Pure functions used by :mod:`epyr.baseline.correction` and the automatic
model-selection code; no data processing logic lives here. All models are
parameterized so that scipy.optimize.curve_fit can call them directly.
"""

from typing import Callable, List, Tuple

import numpy as np


def polynomial_2d(xy: Tuple[np.ndarray, np.ndarray], *coeffs: float) -> np.ndarray:
    """Evaluate a 2D polynomial surface on flattened coordinates.

    The polynomial order is inferred from the number of coefficients:
    a square layout ``(order+1)**2`` is tried first, falling back to a
    rectangular ``(nx+1)*(ny+1)`` decomposition.

    Parameters
    ----------
    xy : tuple of (np.ndarray, np.ndarray)
        Flattened x and y coordinate arrays of identical shape.
    *coeffs : float
        Polynomial coefficients, ordered as ``c[i*(order_y+1)+j]`` for
        the ``x**i * y**j`` term.

    Returns
    -------
    np.ndarray
        Flattened surface values, same shape as ``xy[0]``.

    Raises
    ------
    ValueError
        If ``len(coeffs)`` does not factor as ``(nx+1)*(ny+1)`` with
        ``nx, ny < 10``.

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.baseline.models import polynomial_2d
    >>> x, y = np.meshgrid(np.linspace(0, 1, 5), np.linspace(0, 1, 5))
    >>> z = polynomial_2d((x.ravel(), y.ravel()), 1.0, 2.0, 3.0, 4.0)  # bilinear
    >>> z.shape
    (25,)
    """
    x, y = xy
    result = np.zeros_like(x)
    num_coeffs = len(coeffs)

    # Infer polynomial order from coefficient count
    order = int(np.sqrt(num_coeffs)) - 1
    if (order + 1) ** 2 != num_coeffs:
        # Not square -> search for rectangular factorisation
        for nx in range(10):
            for ny in range(10):
                if (nx + 1) * (ny + 1) == num_coeffs:
                    order_x, order_y = nx, ny
                    break
            else:
                continue
            break
        else:
            raise ValueError(
                f"Cannot determine polynomial order from {num_coeffs} coefficients"
            )
    else:
        order_x = order_y = order

    coeff_idx = 0
    for i in range(order_x + 1):
        for j in range(order_y + 1):
            result += coeffs[coeff_idx] * (x**i) * (y**j)
            coeff_idx += 1

    return result


def stretched_exponential_1d(
    x: np.ndarray,
    A: float,
    tau: float,
    beta: float,
    offset: float = 0.0,
) -> np.ndarray:
    r"""Stretched-exponential (Kohlrausch-Williams-Watts) decay.

    .. math::

        y(x) = \mathrm{offset} + A \exp\!\left[-(x/\tau)^{\beta}\right]

    Commonly used for T2 echo decay in disordered systems.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (time, often in microseconds).
    A : float
        Amplitude.
    tau : float
        Characteristic decay time, same unit as ``x``.
    beta : float
        Stretching exponent in (0, 5]. ``beta = 1`` recovers a pure
        exponential; ``beta < 1`` is sub-exponential; ``beta > 1`` super.
    offset : float, optional
        Constant baseline. Default 0.

    Returns
    -------
    np.ndarray
        Decay values, same shape as ``x``.

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.baseline.models import stretched_exponential_1d
    >>> t = np.linspace(0, 5, 100)  # microseconds
    >>> y = stretched_exponential_1d(t, A=1.0, tau=1.5, beta=0.7)
    >>> y[0], round(y[-1], 4)
    (1.0, 0.098)
    """
    return offset + A * np.exp(-((x / tau) ** beta))


def bi_exponential_1d(
    x: np.ndarray,
    A1: float,
    tau1: float,
    A2: float,
    tau2: float,
    offset: float = 0.0,
) -> np.ndarray:
    r"""Sum of two exponential decays.

    .. math::

        y(x) = \mathrm{offset} + A_1 e^{-x/\tau_1} + A_2 e^{-x/\tau_2}

    Useful when a system has two well-separated relaxation channels.

    Parameters
    ----------
    x : np.ndarray
        Independent variable.
    A1, A2 : float
        Component amplitudes.
    tau1, tau2 : float
        Component decay times, same unit as ``x``. Order does not matter.
    offset : float, optional
        Constant baseline. Default 0.

    Returns
    -------
    np.ndarray
        Decay values.

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.baseline.models import bi_exponential_1d
    >>> t = np.linspace(0, 10, 200)
    >>> y = bi_exponential_1d(t, A1=0.7, tau1=0.5, A2=0.3, tau2=4.0, offset=0.05)
    >>> round(y[0], 3), round(y[-1], 4)
    (1.05, 0.0747)
    """
    return offset + A1 * np.exp(-x / tau1) + A2 * np.exp(-x / tau2)


def polynomial_1d(x: np.ndarray, *coeffs: float) -> np.ndarray:
    r"""Evaluate a 1D polynomial.

    .. math::

        y(x) = a_0 + a_1 x + a_2 x^2 + \ldots + a_n x^n

    Parameters
    ----------
    x : np.ndarray
        Independent variable.
    *coeffs : float
        Coefficients in ascending order of degree, ``[a0, a1, ..., an]``.

    Returns
    -------
    np.ndarray
        Polynomial values (float).

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.baseline.models import polynomial_1d
    >>> x = np.linspace(-1, 1, 5)
    >>> polynomial_1d(x, 1.0, 0.0, 2.0)  # 1 + 2*x**2
    array([3. , 1.5, 1. , 1.5, 3. ])
    """
    result = np.zeros_like(x, dtype=float)
    for i, coeff in enumerate(coeffs):
        result += coeff * (x**i)
    return result


def exponential_1d(
    x: np.ndarray, A: float, tau: float, offset: float = 0.0
) -> np.ndarray:
    r"""Single exponential decay.

    .. math::

        y(x) = \mathrm{offset} + A \exp(-x/\tau)

    Parameters
    ----------
    x : np.ndarray
        Independent variable.
    A : float
        Amplitude.
    tau : float
        Decay time, same unit as ``x``.
    offset : float, optional
        Constant baseline. Default 0.

    Returns
    -------
    np.ndarray
        Decay values.

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.baseline.models import exponential_1d
    >>> t = np.linspace(0, 5, 50)
    >>> y = exponential_1d(t, A=1.0, tau=1.0)
    >>> round(y[0], 3), round(y[-1], 4)
    (1.0, 0.0067)
    """
    return offset + A * np.exp(-x / tau)


# Model information for automatic selection
MODEL_INFO = {
    "polynomial": {
        "name": "Polynomial",
        "function_1d": polynomial_1d,
        "function_2d": polynomial_2d,
        "description": "Polynomial baseline for smooth drifts",
        "typical_use": "CW EPR spectra with baseline drift",
    },
    "exponential": {
        "name": "Exponential",
        "function_1d": exponential_1d,
        "description": "Simple exponential decay",
        "typical_use": "Simple decay processes",
    },
    "stretched_exponential": {
        "name": "Stretched Exponential",
        "function_1d": stretched_exponential_1d,
        "description": "Stretched exponential decay (KWW function)",
        "typical_use": "T2 relaxation, echo decay measurements",
    },
    "bi_exponential": {
        "name": "Bi-exponential",
        "function_1d": bi_exponential_1d,
        "description": "Sum of two exponential decays",
        "typical_use": "Complex decay with multiple components",
    },
}


def get_model_function(
    model_name: str, dimension: str = "1d"
) -> Callable[..., np.ndarray]:
    """Return the callable for a registered model.

    Parameters
    ----------
    model_name : str
        Key from :data:`MODEL_INFO` (``'polynomial'``,
        ``'stretched_exponential'``, ``'bi_exponential'``, ``'exponential'``).
    dimension : {'1d', '2d'}, optional
        Function variant. Most models only provide ``'1d'``; only the
        polynomial model supports ``'2d'``. Default ``'1d'``.

    Returns
    -------
    callable
        The model function, ready to pass to ``scipy.optimize.curve_fit``.

    Raises
    ------
    ValueError
        If ``model_name`` is unknown or the requested ``dimension`` variant
        does not exist.

    Examples
    --------
    >>> from epyr.baseline.models import get_model_function
    >>> f = get_model_function("stretched_exponential")
    >>> f.__name__
    'stretched_exponential_1d'
    """
    if model_name not in MODEL_INFO:
        raise ValueError(f"Unknown model: {model_name}")

    func_key = f"function_{dimension}"
    if func_key not in MODEL_INFO[model_name]:
        raise ValueError(f"Model {model_name} does not support {dimension}")

    return MODEL_INFO[model_name][func_key]


def list_available_models() -> List[str]:
    """Return the names of registered baseline models.

    Returns
    -------
    list of str
        Model keys usable with :func:`get_model_function` and
        :data:`MODEL_INFO`.

    Examples
    --------
    >>> from epyr.baseline.models import list_available_models
    >>> sorted(list_available_models())
    ['bi_exponential', 'exponential', 'polynomial', 'stretched_exponential']
    """
    return list(MODEL_INFO.keys())
