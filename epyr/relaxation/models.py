"""
Time-domain relaxation models for EPR T1/T2 measurements.

Mathematical Background:
    Mono-exponential decay/recovery, stretched exponential, bi-exponential,
    and a combined homogeneous/spectral-diffusion decay model used to
    extract spin relaxation times (T1, T2) from pulsed EPR time traces.

References:
    Eaton, S.S. and Eaton, G.R., "Relaxation Times of Organic Radicals and
    Transition Metal Ions", Biol. Magn. Reson., 19, 2000.
"""

import numpy as np


def mono_exponential(
    t: np.ndarray,
    amplitude: float,
    T: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a single-rate exponential decay.

    V(t) = amplitude * exp(-t / T) + offset

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller (e.g. ns or us).
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    T : float
        Decay (or recovery) time constant, in the same unit as t.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * np.exp(-t / T) + offset


def stretched_exponential(
    t: np.ndarray,
    amplitude: float,
    T: float,
    beta: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a stretched-exponential decay.

    V(t) = amplitude * exp(-(t / T) ** beta) + offset

    beta = 1 reduces to a mono-exponential decay. beta < 1 describes a
    distribution of relaxation rates, common in pulsed EPR echo decays.

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    T : float
        Characteristic decay time constant, in the same unit as t.
    beta : float
        Stretching exponent (dimensionless).
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * np.exp(-((t / T) ** beta)) + offset


def biexponential(
    t: np.ndarray,
    amplitude1: float,
    tau1: float,
    amplitude2: float,
    tau2: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a two-component exponential decay.

    V(t) = amplitude1 * exp(-t / tau1) + amplitude2 * exp(-t / tau2) + offset

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude1 : float
        Amplitude of the first component at t = 0.
    tau1 : float
        Decay time constant of the first component, same unit as t.
    amplitude2 : float
        Amplitude of the second component at t = 0.
    tau2 : float
        Decay time constant of the second component, same unit as t.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.

    Notes
    -----
    tau1 and tau2 are not necessarily T1 or T2: either component can
    represent either relaxation process, depending on the experiment.
    """
    t = np.asarray(t, dtype=float)
    return (
        amplitude1 * np.exp(-t / tau1) + amplitude2 * np.exp(-t / tau2) + offset
    )


def inversion_recovery(
    t: np.ndarray,
    amplitude: float,
    T1: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute an inversion-recovery T1 curve.

    V(t) = amplitude * (1 - 2 * exp(-t / T1)) + offset

    Parameters
    ----------
    t : np.ndarray
        Delay time after inversion, in the unit chosen by the caller.
    amplitude : float
        Fully-recovered signal amplitude relative to offset.
    T1 : float
        Spin-lattice relaxation time, same unit as t.
    offset : float, optional
        Baseline signal level (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * (1.0 - 2.0 * np.exp(-t / T1)) + offset


def saturation_recovery(
    t: np.ndarray,
    amplitude: float,
    T1: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a saturation-recovery T1 curve.

    V(t) = amplitude * (1 - exp(-t / T1)) + offset

    Parameters
    ----------
    t : np.ndarray
        Delay time after saturation, in the unit chosen by the caller.
    amplitude : float
        Fully-recovered signal amplitude relative to offset.
    T1 : float
        Spin-lattice relaxation time, same unit as t.
    offset : float, optional
        Baseline signal level (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * (1.0 - np.exp(-t / T1)) + offset


def gamma_gaussian_decay(
    t: np.ndarray,
    amplitude: float,
    Gamma0: float,
    GammaG: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute an echo decay with homogeneous and spectral-diffusion terms.

    V(t) = amplitude * exp(-Gamma0 * t) * exp(-(GammaG * t) ** 2) + offset

    Gamma0 captures the homogeneous (exponential) relaxation rate. GammaG
    captures an additional Gaussian-shaped decay from spectral diffusion.

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    Gamma0 : float
        Homogeneous decay rate, in inverse time units (1 / t unit).
    GammaG : float
        Gaussian (spectral-diffusion) decay rate, in inverse time units.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return (
        amplitude * np.exp(-Gamma0 * t) * np.exp(-((GammaG * t) ** 2)) + offset
    )
