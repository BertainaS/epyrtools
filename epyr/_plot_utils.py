"""
Shared matplotlib helpers for the fit-plotting code in
:mod:`epyr.lineshapes.fitting` and :mod:`epyr.relaxation.fitting`.
"""

from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def get_or_create_subplots(
    nrows: int,
    ncols: int,
    gridspec_kw: Optional[Dict[str, Any]] = None,
) -> Tuple[Figure, Any]:
    """
    Return subplots on the current figure if one is open, else a new figure.

    Reusing an open figure lets repeated calls (e.g. re-running a fit in a
    notebook cell) update the same window in place instead of accumulating
    one figure per call. Figure size and other style defaults come from
    ``matplotlib.rcParams`` either way.

    Parameters
    ----------
    nrows : int
        Number of subplot rows.
    ncols : int
        Number of subplot columns.
    gridspec_kw : dict, optional
        Forwarded to ``Figure.subplots`` / ``pyplot.subplots``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The reused or newly created figure.
    axes : matplotlib.axes.Axes or array of Axes
        Subplot axes, in the same layout ``pyplot.subplots`` would return.
    """
    if plt.get_fignums():
        fig = plt.gcf()
        fig.clf()
        axes = fig.subplots(nrows, ncols, gridspec_kw=gridspec_kw)
    else:
        fig, axes = plt.subplots(nrows, ncols, gridspec_kw=gridspec_kw)

    return fig, axes
