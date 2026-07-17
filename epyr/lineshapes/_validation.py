"""Shared input validation for lineshape functions."""

import numpy as np


def validate_abscissa(x) -> np.ndarray:
    """Convert abscissa input to a float array and reject empty input.

    Parameters
    ----------
    x : array_like
        Abscissa points, in the same unit as the lineshape width.

    Returns
    -------
    np.ndarray
        x as a float64 array.

    Raises
    ------
    ValueError
        If x contains no points.
    """
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        raise ValueError("x must contain at least one point")
    return x
