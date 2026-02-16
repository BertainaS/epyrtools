#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Region selection utilities for baseline correction.

This module provides functions for creating and managing baseline regions
for both 1D and 2D EPR data.
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


def create_region_mask_1d(
    x: np.ndarray, regions: List[Tuple[float, float]], mode: str = "exclude"
) -> np.ndarray:
    """
    Create a boolean mask for 1D data based on selected regions.

    Args:
        x: X-coordinate array
        regions: List of regions as [(x1, x2), ...]
        mode: 'exclude' to exclude regions, 'include' to include only regions

    Returns:
        Boolean mask array (True = use for fitting)
    """
    if mode == "exclude":
        mask = np.ones(len(x), dtype=bool)
        for x1, x2 in regions:
            mask &= ~((x >= min(x1, x2)) & (x <= max(x1, x2)))
    elif mode == "include":
        mask = np.zeros(len(x), dtype=bool)
        for x1, x2 in regions:
            mask |= (x >= min(x1, x2)) & (x <= max(x1, x2))
    else:
        raise ValueError("mode must be 'exclude' or 'include'")

    return mask


def create_region_mask_2d(
    X: np.ndarray,
    Y: np.ndarray,
    regions: List[Tuple[Tuple[float, float], Tuple[float, float]]],
    mode: str = "exclude",
) -> np.ndarray:
    """
    Create a boolean mask for 2D data based on selected regions.

    Args:
        X: X-coordinate meshgrid
        Y: Y-coordinate meshgrid
        regions: List of regions as [((x1,x2), (y1,y2)), ...]
        mode: 'exclude' to exclude regions, 'include' to include only regions

    Returns:
        Boolean mask array (True = use for fitting)
    """
    if mode == "exclude":
        mask = np.ones_like(X, dtype=bool)
        for (x1, x2), (y1, y2) in regions:
            x1, x2 = min(x1, x2), max(x1, x2)
            y1, y2 = min(y1, y2), max(y1, y2)
            mask &= ~((X >= x1) & (X <= x2) & (Y >= y1) & (Y <= y2))
    elif mode == "include":
        mask = np.zeros_like(X, dtype=bool)
        for (x1, x2), (y1, y2) in regions:
            x1, x2 = min(x1, x2), max(x1, x2)
            y1, y2 = min(y1, y2), max(y1, y2)
            mask |= (X >= x1) & (X <= x2) & (Y >= y1) & (Y <= y2)
    else:
        raise ValueError("mode must be 'exclude' or 'include'")

    return mask


def create_center_exclusion_mask_1d(
    x: np.ndarray, center_fraction: float = 0.3
) -> np.ndarray:
    """
    Create a mask that excludes the center portion of 1D data.

    This is useful for CW EPR spectra where the signal is in the center
    and we want to fit the baseline using the wings.

    Args:
        x: X-coordinate array
        center_fraction: Fraction of the data range to exclude from center

    Returns:
        Boolean mask array (True = use for fitting, False = exclude)
    """
    x_min, x_max = x.min(), x.max()
    x_range = x_max - x_min
    center = (x_min + x_max) / 2

    exclude_half_width = (center_fraction * x_range) / 2
    exclude_min = center - exclude_half_width
    exclude_max = center + exclude_half_width

    mask = (x < exclude_min) | (x > exclude_max)
    return mask


def create_center_exclusion_mask_2d(
    X: np.ndarray, Y: np.ndarray, center_fraction: float = 0.3
) -> np.ndarray:
    """
    Create a mask that excludes the center portion of 2D data.

    Args:
        X: X-coordinate meshgrid
        Y: Y-coordinate meshgrid
        center_fraction: Fraction of the data range to exclude from center

    Returns:
        Boolean mask array (True = use for fitting, False = exclude)
    """
    # X-direction exclusion
    x_min, x_max = X.min(), X.max()
    x_range = x_max - x_min
    x_center = (x_min + x_max) / 2
    x_exclude_half_width = (center_fraction * x_range) / 2

    # Y-direction exclusion
    y_min, y_max = Y.min(), Y.max()
    y_range = y_max - y_min
    y_center = (y_min + y_max) / 2
    y_exclude_half_width = (center_fraction * y_range) / 2

    # Create exclusion mask (exclude center rectangle)
    x_mask = (X < (x_center - x_exclude_half_width)) | (
        X > (x_center + x_exclude_half_width)
    )
    y_mask = (Y < (y_center - y_exclude_half_width)) | (
        Y > (y_center + y_exclude_half_width)
    )

    # Include point if it's outside the center region in either x OR y
    mask = x_mask | y_mask
    return mask


def create_edge_exclusion_mask_1d(
    x: np.ndarray, exclude_initial: int = 0, exclude_final: int = 0
) -> np.ndarray:
    """
    Create a mask that excludes points at the beginning and end of 1D data.

    This is useful for time-series data where initial/final points may
    have artifacts or noise.

    Args:
        x: X-coordinate array
        exclude_initial: Number of initial points to exclude
        exclude_final: Number of final points to exclude

    Returns:
        Boolean mask array (True = use for fitting)
    """
    mask = np.ones(len(x), dtype=bool)

    if exclude_initial > 0:
        mask[:exclude_initial] = False
    if exclude_final > 0:
        mask[-exclude_final:] = False

    return mask


def combine_masks(*masks: np.ndarray) -> np.ndarray:
    """
    Combine multiple boolean masks using logical AND.

    Args:
        masks: Variable number of boolean mask arrays

    Returns:
        Combined boolean mask (True where ALL masks are True)
    """
    if not masks:
        raise ValueError("At least one mask must be provided")

    combined = masks[0].copy()
    for mask in masks[1:]:
        combined &= mask

    return combined


def get_baseline_regions_1d(
    x: np.ndarray,
    y: np.ndarray,
    exclude_center: bool = True,
    center_fraction: float = 0.3,
    exclude_initial: int = 0,
    exclude_final: int = 0,
    manual_regions: Optional[List[Tuple[float, float]]] = None,
    region_mode: str = "exclude",
) -> np.ndarray:
    """
    Create a comprehensive baseline mask for 1D data combining all exclusion criteria.

    Args:
        x: X-coordinate array
        y: Y-data array
        exclude_center: Whether to exclude center region
        center_fraction: Fraction of center to exclude
        exclude_initial: Number of initial points to exclude
        exclude_final: Number of final points to exclude
        manual_regions: List of manually specified regions
        region_mode: How to handle manual_regions ('exclude' or 'include')

    Returns:
        Boolean mask array for baseline fitting
    """
    masks = []

    # Edge exclusion mask
    if exclude_initial > 0 or exclude_final > 0:
        edge_mask = create_edge_exclusion_mask_1d(x, exclude_initial, exclude_final)
        masks.append(edge_mask)

    # Center exclusion mask
    if exclude_center:
        center_mask = create_center_exclusion_mask_1d(x, center_fraction)
        masks.append(center_mask)

    # Manual regions mask
    if manual_regions:
        manual_mask = create_region_mask_1d(x, manual_regions, region_mode)
        masks.append(manual_mask)

    # Combine all masks
    if masks:
        return combine_masks(*masks)
    else:
        # No exclusions, use all points
        return np.ones(len(x), dtype=bool)


def get_baseline_regions_2d(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    exclude_center: bool = True,
    center_fraction: float = 0.3,
    manual_regions: Optional[
        List[Tuple[Tuple[float, float], Tuple[float, float]]]
    ] = None,
    region_mode: str = "exclude",
) -> np.ndarray:
    """
    Create a comprehensive baseline mask for 2D data combining all exclusion criteria.

    Args:
        X: X-coordinate meshgrid
        Y: Y-coordinate meshgrid
        Z: Z-data array
        exclude_center: Whether to exclude center region
        center_fraction: Fraction of center to exclude
        manual_regions: List of manually specified regions
        region_mode: How to handle manual_regions ('exclude' or 'include')

    Returns:
        Boolean mask array for baseline fitting
    """
    masks = []

    # Center exclusion mask
    if exclude_center:
        center_mask = create_center_exclusion_mask_2d(X, Y, center_fraction)
        masks.append(center_mask)

    # Manual regions mask
    if manual_regions:
        manual_mask = create_region_mask_2d(X, Y, manual_regions, region_mode)
        masks.append(manual_mask)

    # Combine all masks
    if masks:
        return combine_masks(*masks)
    else:
        # No exclusions, use all points
        return np.ones_like(X, dtype=bool)


def validate_regions_1d(
    regions: List[Tuple[float, float]], x_min: float, x_max: float
) -> bool:
    """
    Validate that 1D regions are within data bounds and properly formatted.

    Args:
        regions: List of region tuples
        x_min: Minimum x value in data
        x_max: Maximum x value in data

    Returns:
        True if all regions are valid

    Raises:
        ValueError: If regions are invalid
    """
    for i, (x1, x2) in enumerate(regions):
        if not (x_min <= x1 <= x_max and x_min <= x2 <= x_max):
            raise ValueError(
                f"Region {i} ({x1}, {x2}) is outside data bounds ({x_min}, {x_max})"
            )
        if x1 == x2:
            raise ValueError(f"Region {i} has zero width: ({x1}, {x2})")

    return True


def validate_regions_2d(
    regions: List[Tuple[Tuple[float, float], Tuple[float, float]]],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> bool:
    """
    Validate that 2D regions are within data bounds and properly formatted.

    Args:
        regions: List of region tuples
        x_min, x_max: X data bounds
        y_min, y_max: Y data bounds

    Returns:
        True if all regions are valid

    Raises:
        ValueError: If regions are invalid
    """
    for i, ((x1, x2), (y1, y2)) in enumerate(regions):
        if not (x_min <= x1 <= x_max and x_min <= x2 <= x_max):
            raise ValueError(
                f"Region {i} X bounds ({x1}, {x2}) outside data bounds ({x_min}, {x_max})"
            )
        if not (y_min <= y1 <= y_max and y_min <= y2 <= y_max):
            raise ValueError(
                f"Region {i} Y bounds ({y1}, {y2}) outside data bounds ({y_min}, {y_max})"
            )
        if x1 == x2 or y1 == y2:
            raise ValueError(f"Region {i} has zero area: (({x1}, {x2}), ({y1}, {y2}))")

    return True
