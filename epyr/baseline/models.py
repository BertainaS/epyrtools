#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mathematical models for baseline correction.

This module contains the mathematical functions used for baseline fitting
in EPR data analysis. These are pure mathematical functions without any
data processing logic.
"""

from typing import Tuple

import numpy as np


def polynomial_2d(xy, *coeffs):
    """
    2D polynomial function for curve_fit.

    Args:
        xy: tuple of (x, y) coordinate arrays (flattened)
        coeffs: polynomial coefficients

    Returns:
        Flattened polynomial surface values
    """
    x, y = xy
    result = np.zeros_like(x)

    # Determine polynomial order from number of coefficients
    # For order (nx, ny): num_coeffs = (nx+1) * (ny+1)
    num_coeffs = len(coeffs)

    # Find polynomial orders (assume square for simplicity, or deduce from coeffs)
    # For now, assume same order in both directions
    order = int(np.sqrt(num_coeffs)) - 1
    if (order + 1) ** 2 != num_coeffs:
        # If not square, try rectangular
        for nx in range(10):  # reasonable limit
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

    # Calculate polynomial
    coeff_idx = 0
    for i in range(order_x + 1):
        for j in range(order_y + 1):
            result += coeffs[coeff_idx] * (x**i) * (y**j)
            coeff_idx += 1

    return result


def stretched_exponential_1d(x, A, tau, beta, offset=0):
    """
    Stretched exponential decay function for baseline fitting.

    This function represents a stretched exponential decay commonly observed
    in EPR relaxation measurements, particularly T2 echo decay data.

    Args:
        x: Time or field array
        A: Amplitude (positive for decay)
        tau: Decay time constant (characteristic time scale)
        beta: Stretching exponent (0 < beta <= 5.0)
              - beta = 1: Pure exponential decay
              - beta < 1: Sub-exponential (slower than exponential)
              - beta > 1: Super-exponential (faster than exponential)
        offset: Constant offset

    Returns:
        np.ndarray: Stretched exponential values = offset + A * exp(-(x/tau)^beta)
    """
    return offset + A * np.exp(-((x / tau) ** beta))


def bi_exponential_1d(x, A1, tau1, A2, tau2, offset=0):
    """
    Bi-exponential decay function for baseline fitting.

    This function represents the sum of two exponential decays, useful for
    systems with multiple relaxation pathways or decay components.

    Args:
        x: Time or field array
        A1: Amplitude of first exponential component
        tau1: Decay time constant of first exponential
        A2: Amplitude of second exponential component
        tau2: Decay time constant of second exponential
        offset: Constant offset

    Returns:
        np.ndarray: Bi-exponential values = offset + A1*exp(-x/tau1) + A2*exp(-x/tau2)
    """
    return offset + A1 * np.exp(-x / tau1) + A2 * np.exp(-x / tau2)


def polynomial_1d(x, *coeffs):
    """
    1D polynomial function for baseline fitting.

    Args:
        x: Independent variable array
        coeffs: Polynomial coefficients [a0, a1, a2, ..., an]
                Result = a0 + a1*x + a2*x^2 + ... + an*x^n

    Returns:
        np.ndarray: Polynomial values
    """
    result = np.zeros_like(x)
    for i, coeff in enumerate(coeffs):
        result += coeff * (x**i)
    return result


def exponential_1d(x, A, tau, offset=0):
    """
    Simple exponential decay function.

    Args:
        x: Time or field array
        A: Amplitude
        tau: Decay time constant
        offset: Constant offset

    Returns:
        np.ndarray: Exponential values = offset + A * exp(-x/tau)
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


def get_model_function(model_name: str, dimension: str = "1d"):
    """
    Get the mathematical function for a given model.

    Args:
        model_name: Name of the model ('polynomial', 'stretched_exponential', etc.)
        dimension: '1d' or '2d'

    Returns:
        Callable: Mathematical function
    """
    if model_name not in MODEL_INFO:
        raise ValueError(f"Unknown model: {model_name}")

    func_key = f"function_{dimension}"
    if func_key not in MODEL_INFO[model_name]:
        raise ValueError(f"Model {model_name} does not support {dimension}")

    return MODEL_INFO[model_name][func_key]


def list_available_models() -> list:
    """Get list of available baseline models."""
    return list(MODEL_INFO.keys())
