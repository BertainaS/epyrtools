#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automatic baseline model selection for EPR data.

This module provides intelligent automatic selection of the best baseline
correction model using statistical criteria (AIC, BIC, R²).
"""

import sys
import warnings
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..logging_config import get_logger

logger = get_logger(__name__)

# Import correction functions from our new modules
from .correction import (
    baseline_bi_exponential_1d,
    baseline_polynomial_1d,
    baseline_stretched_exponential_1d,
)
from .interactive import RegionSelector, is_interactive_available
from .models import MODEL_INFO


def _calculate_model_metrics(y_data, y_pred, n_params, n_points):
    """
    Calculate model selection metrics (AIC, BIC, R²).

    Args:
        y_data: Observed data
        y_pred: Model predictions
        n_params: Number of model parameters
        n_points: Number of data points

    Returns:
        dict: Dictionary with 'aic', 'bic', 'r2', 'rss'
    """
    residuals = y_data - y_pred
    rss = np.sum(residuals**2)

    # Avoid division by zero
    if rss == 0:
        rss = 1e-10

    # R-squared
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r2 = 1 - rss / ss_tot if ss_tot > 0 else 0

    # AIC and BIC
    try:
        aic = 2 * n_params + n_points * np.log(rss / n_points)
        bic = n_params * np.log(n_points) + n_points * np.log(rss / n_points)
    except (ValueError, OverflowError):
        # Handle edge cases
        aic = np.inf
        bic = np.inf

    return {"aic": aic, "bic": bic, "r2": r2, "rss": rss}


def _test_polynomial_models(
    x, y, params, common_args, selection_criterion="aic", verbose=True
):
    """Test polynomial models of different orders and return the best one."""
    best_result = None
    best_criterion = np.inf if selection_criterion in ["aic", "bic"] else -np.inf

    for order in [1, 2, 3, 4]:
        try:
            # Suppress output for cleaner comparison
            old_stdout = sys.stdout
            sys.stdout = StringIO()

            try:
                # Polynomial functions have different parameter names
                poly_args = {
                    "manual_regions": common_args.get("manual_regions"),
                    "region_mode": common_args.get("region_mode"),
                    "interactive": False,
                    "exclude_center": (
                        True if not common_args.get("manual_regions") else False
                    ),
                    "center_fraction": 0.3,
                }
                # Remove None values
                poly_args = {k: v for k, v in poly_args.items() if v is not None}

                corrected, baseline = baseline_polynomial_1d(
                    x, y, params, order=order, **poly_args
                )

            finally:
                sys.stdout = old_stdout

            # Calculate metrics
            metrics = _calculate_model_metrics(y, baseline, order + 1, len(y))
            criterion_value = metrics[selection_criterion]

            # Check if this is the best so far
            is_better = (
                selection_criterion in ["aic", "bic"]
                and criterion_value < best_criterion
            ) or (selection_criterion == "r2" and criterion_value > best_criterion)

            if is_better:
                best_criterion = criterion_value
                best_result = {
                    "corrected": corrected,
                    "baseline": baseline,
                    "order": order,
                    "n_params": order + 1,
                    "metrics": metrics,
                }

        except Exception as e:
            if verbose:
                logger.warning(f"     Polynomial order {order}: FAILED - {e}")
            continue

    return best_result


def _test_stretched_exponential_model(
    x, y, params, common_args, selection_criterion="aic", verbose=True
):
    """Test stretched exponential model."""
    try:
        # Suppress output for cleaner comparison
        old_stdout = sys.stdout
        sys.stdout = StringIO()

        try:
            corrected, baseline = baseline_stretched_exponential_1d(
                x, y, params, **common_args
            )
        finally:
            sys.stdout = old_stdout

        # Calculate metrics (4 parameters: A, tau, beta, offset)
        metrics = _calculate_model_metrics(y, baseline, 4, len(y))

        result = {
            "corrected": corrected,
            "baseline": baseline,
            "n_params": 4,
            "metrics": metrics,
        }

        return result

    except Exception as e:
        if verbose:
            logger.warning(f"   Stretched exponential model failed: {e}")
        return None


def _test_bi_exponential_model(
    x, y, params, common_args, selection_criterion="aic", verbose=True
):
    """Test bi-exponential model."""
    try:
        # Suppress output for cleaner comparison
        old_stdout = sys.stdout
        sys.stdout = StringIO()

        try:
            corrected, baseline = baseline_bi_exponential_1d(
                x, y, params, **common_args
            )
        finally:
            sys.stdout = old_stdout

        # Calculate metrics (5 parameters: A1, tau1, A2, tau2, offset)
        metrics = _calculate_model_metrics(y, baseline, 5, len(y))

        result = {
            "corrected": corrected,
            "baseline": baseline,
            "n_params": 5,
            "metrics": metrics,
        }

        return result

    except Exception as e:
        if verbose:
            logger.warning(f"   Bi-exponential model failed: {e}")
        return None


def baseline_auto_1d(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    models: List[str] = ["polynomial", "stretched_exponential", "bi_exponential"],
    selection_criterion: str = "aic",
    use_real_part: bool = True,
    exclude_initial: int = 0,
    exclude_final: int = 0,
    manual_regions: Optional[List[Tuple[float, float]]] = None,
    region_mode: str = "exclude",
    interactive: bool = False,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Automatic baseline model selection and correction for 1D EPR data.

    This function tests multiple baseline models and automatically selects
    the best one based on information criteria (AIC/BIC) or other metrics.

    Args:
        x: X-axis data from eprload (can be None)
        y: 1D spectral data array from eprload
        params: Parameter dictionary from eprload (optional)
        models: List of models to test ['polynomial', 'stretched_exponential', 'bi_exponential']
        selection_criterion: 'aic' (Akaike), 'bic' (Bayesian), or 'r2' (R-squared)
        use_real_part: If True, fit only real part of complex data
        exclude_initial: Number of initial points to exclude from fitting
        exclude_final: Number of final points to exclude from fitting
        manual_regions: List of manually selected regions as [(x1, x2), ...]
        region_mode: 'exclude' to exclude manual_regions, 'include' to use only manual_regions
        interactive: If True, open interactive region selector
        verbose: If True, print detailed model comparison

    Returns:
        tuple: (corrected_data, best_baseline, model_info)
            model_info contains: {'best_model': str, 'criteria': dict, 'parameters': dict}

    Examples:
        # Automatic model selection
        corrected, baseline, info = baseline_auto_1d(x, y, params)
        print(f"Best model: {info['best_model']}")

        # Restrict to specific models
        corrected, baseline, info = baseline_auto_1d(x, y, params, models=['polynomial', 'stretched_exponential'])

        # Use BIC for model selection
        corrected, baseline, info = baseline_auto_1d(x, y, params, selection_criterion='bic')
    """
    if y is None or y.ndim != 1:
        raise ValueError("y must be a 1D array")

    if selection_criterion not in ["aic", "bic", "r2"]:
        raise ValueError("selection_criterion must be 'aic', 'bic', or 'r2'")

    n_points = len(y)

    # Handle complex data
    if np.iscomplexobj(y):
        if use_real_part:
            y_fit = np.real(y)
            if verbose:
                logger.info("ℹ Using real part of complex data for fitting")
        else:
            y_fit = np.abs(y)
            if verbose:
                logger.info("ℹ Using magnitude of complex data for fitting")
    else:
        y_fit = y.copy()

    # Create x-coordinates
    if x is None or len(x) != n_points:
        x_coords = np.arange(n_points, dtype=float)
        if verbose:
            logger.info("ℹ No x-axis provided, using index as x-coordinates")
    else:
        x_coords = x.astype(float)

    # Interactive region selection if requested
    selected_regions = manual_regions if manual_regions is not None else []

    if interactive:
        if not is_interactive_available():
            logger.warning("⚠️  Interactive selection may not work in this environment.")

        if verbose:
            logger.info("🖱️ Interactive region selection enabled...")

        selector = RegionSelector()
        selected_regions = selector.select_regions_1d(
            x_coords,
            y_fit,
            f"Select regions to {region_mode.upper()} from automatic baseline fitting",
        )
        if verbose:
            logger.info(f"✅ Selected {len(selected_regions)} regions")

    # Prepare common arguments for all models
    common_args = {
        "use_real_part": use_real_part,
        "exclude_initial": exclude_initial,
        "exclude_final": exclude_final,
        "manual_regions": selected_regions,
        "region_mode": region_mode,
        "interactive": False,  # We already did interactive selection above
    }

    # Test each model
    model_results = {}

    if verbose:
        logger.info(f"\n🧪 Testing {len(models)} baseline models...")

    # Test polynomial model
    if "polynomial" in models:
        if verbose:
            logger.info("   Testing polynomial baseline...")

        poly_result = _test_polynomial_models(
            x_coords, y_fit, params, common_args, selection_criterion, verbose
        )
        if poly_result is not None:
            model_results["polynomial"] = poly_result
            if verbose:
                metrics = poly_result["metrics"]
                logger.info(
                    f"   ✅ Polynomial (order {poly_result['order']}): {selection_criterion.upper()}={metrics[selection_criterion]:.2f}"
                )

    # Test stretched exponential model
    if "stretched_exponential" in models:
        if verbose:
            logger.info("   Testing stretched exponential baseline...")

        stretch_result = _test_stretched_exponential_model(
            x_coords, y_fit, params, common_args, selection_criterion, verbose
        )
        if stretch_result is not None:
            model_results["stretched_exponential"] = stretch_result
            if verbose:
                metrics = stretch_result["metrics"]
                logger.info(
                    f"   ✅ Stretched exponential: {selection_criterion.upper()}={metrics[selection_criterion]:.2f}"
                )

    # Test bi-exponential model
    if "bi_exponential" in models:
        if verbose:
            logger.info("   Testing bi-exponential baseline...")

        bi_result = _test_bi_exponential_model(
            x_coords, y_fit, params, common_args, selection_criterion, verbose
        )
        if bi_result is not None:
            model_results["bi_exponential"] = bi_result
            if verbose:
                metrics = bi_result["metrics"]
                logger.info(
                    f"   ✅ Bi-exponential: {selection_criterion.upper()}={metrics[selection_criterion]:.2f}"
                )

    # Select best model
    if not model_results:
        raise ValueError("All baseline models failed. Check your data and parameters.")

    # Find the best model based on the selection criterion
    best_model = None
    best_criterion = np.inf if selection_criterion in ["aic", "bic"] else -np.inf

    criteria_dict = {}
    for model_name, result in model_results.items():
        criterion_value = result["metrics"][selection_criterion]
        criteria_dict[model_name] = criterion_value

        is_better = (
            selection_criterion in ["aic", "bic"] and criterion_value < best_criterion
        ) or (selection_criterion == "r2" and criterion_value > best_criterion)

        if is_better:
            best_criterion = criterion_value
            best_model = model_name

    if best_model is None:
        raise ValueError("Could not select best baseline model")

    # Get the best result
    best_result = model_results[best_model]

    # Prepare return information
    model_info = {
        "best_model": best_model,
        "criteria": criteria_dict,
        "parameters": best_result["metrics"].copy(),
        "n_params": best_result["n_params"],
    }

    # Add model-specific information
    if best_model == "polynomial":
        model_info["polynomial_order"] = best_result["order"]

    if verbose:
        logger.info(f"\n🏆 Best model: {best_model}")
        logger.info(f"📊 {selection_criterion.upper()} = {best_criterion:.2f}")
        logger.info(f"📊 R² = {best_result['metrics']['r2']:.4f}")

        if len(model_results) > 1:
            logger.info(f"\n📋 Model comparison ({selection_criterion.upper()}):")
            sorted_models = sorted(
                criteria_dict.items(),
                key=lambda x: x[1],
                reverse=(selection_criterion == "r2"),
            )
            for i, (model, criterion_val) in enumerate(sorted_models):
                marker = "🥇" if model == best_model else f"{i+1}. "
                logger.info(f"   {marker} {model}: {criterion_val:.2f}")

    # Return the corrected data in the original format
    corrected_data = best_result["corrected"]
    baseline = best_result["baseline"]

    return corrected_data, baseline, model_info


def compare_models_detailed(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    models: List[str] = ["polynomial", "stretched_exponential", "bi_exponential"],
    **kwargs,
) -> Dict[str, Dict[str, Any]]:
    """
    Detailed comparison of all baseline models with full metrics.

    Args:
        x: X-axis data
        y: Y-data array
        params: Parameter dictionary
        models: List of models to compare
        **kwargs: Additional arguments passed to baseline correction functions

    Returns:
        dict: Detailed results for each model
    """
    # Remove verbose and selection_criterion from kwargs for individual model testing
    kwargs_clean = kwargs.copy()
    kwargs_clean.pop("verbose", None)
    kwargs_clean.pop("selection_criterion", None)

    results = {}

    # Prepare data
    y_fit = y
    if np.iscomplexobj(y):
        if kwargs.get("use_real_part", True):
            y_fit = np.real(y)
        else:
            y_fit = np.abs(y)

    n_points = len(y_fit)
    x_coords = x if x is not None else np.arange(n_points, dtype=float)

    # Common arguments
    common_args = {
        "use_real_part": kwargs.get("use_real_part", True),
        "exclude_initial": kwargs.get("exclude_initial", 0),
        "exclude_final": kwargs.get("exclude_final", 0),
        "manual_regions": kwargs.get("manual_regions"),
        "region_mode": kwargs.get("region_mode", "exclude"),
        "interactive": False,
    }

    # Test each model
    if "polynomial" in models:
        poly_result = _test_polynomial_models(
            x_coords, y_fit, params, common_args, "aic", False
        )
        if poly_result:
            results["polynomial"] = poly_result

    if "stretched_exponential" in models:
        stretch_result = _test_stretched_exponential_model(
            x_coords, y_fit, params, common_args, "aic", False
        )
        if stretch_result:
            results["stretched_exponential"] = stretch_result

    if "bi_exponential" in models:
        bi_result = _test_bi_exponential_model(
            x_coords, y_fit, params, common_args, "aic", False
        )
        if bi_result:
            results["bi_exponential"] = bi_result

    return results


def get_model_recommendations(
    data_type: str = None, experiment_type: str = None
) -> List[str]:
    """
    Get recommended models based on data characteristics.

    Args:
        data_type: 'cw' for continuous wave, 'pulsed' for pulsed EPR
        experiment_type: 't1', 't2', 'echo', 'rabi', etc.

    Returns:
        list: Recommended models in order of preference
    """
    if data_type == "cw":
        return ["polynomial"]

    elif data_type == "pulsed":
        if experiment_type in ["t1", "t2", "echo"]:
            return ["stretched_exponential", "bi_exponential", "polynomial"]
        elif experiment_type in ["rabi", "nutation"]:
            return ["polynomial"]
        else:
            return ["stretched_exponential", "polynomial", "bi_exponential"]

    else:
        # General case - try all models
        return ["polynomial", "stretched_exponential", "bi_exponential"]


def auto_baseline_with_recommendations(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    data_type: str = None,
    experiment_type: str = None,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Automatic baseline correction with intelligent model recommendations.

    Args:
        x: X-axis data
        y: Y-data array
        params: Parameter dictionary
        data_type: Type of EPR data ('cw' or 'pulsed')
        experiment_type: Specific experiment type ('t1', 't2', 'echo', 'rabi', etc.)
        **kwargs: Additional arguments

    Returns:
        tuple: (corrected_data, baseline, model_info)
    """
    recommended_models = get_model_recommendations(data_type, experiment_type)

    kwargs["models"] = recommended_models
    kwargs.setdefault("verbose", True)

    if kwargs.get("verbose"):
        logger.info(
            f"🎯 Recommended models for {data_type or 'unknown'} {experiment_type or 'data'}: {recommended_models}"
        )

    return baseline_auto_1d(x, y, params, **kwargs)
