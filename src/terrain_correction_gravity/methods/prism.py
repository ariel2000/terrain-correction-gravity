"""Exact vertical attraction of right rectangular prisms after Nagy (1966)."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, SI_TO_MGAL, G
from ..io.terrain import TerrainCells


def _safe_xlog(x: np.ndarray, value: np.ndarray) -> np.ndarray:
    result = np.zeros_like(value, dtype=np.float64)
    valid = (x != 0.0) & (value > 0.0)
    result[valid] = x[valid] * np.log(value[valid])
    return result


def vertical_prism_acceleration(
    x1: np.ndarray,
    x2: np.ndarray,
    y1: np.ndarray,
    y2: np.ndarray,
    z1: np.ndarray,
    z2: np.ndarray,
    density_kg_m3: np.ndarray | float,
) -> np.ndarray:
    """Return upward acceleration (m/s^2), with z positive upward.

    The arctangent is evaluated as arctan(x*y/(z*r)), matching the continuous
    Nagy expression on either side of the observation plane. Indeterminate
    zero-limit terms are explicitly set to zero.
    """
    total = np.zeros(np.broadcast(x1, x2, y1, y2, z1, z2).shape, dtype=np.float64)
    for i, x in enumerate((x1, x2)):
        for j, y in enumerate((y1, y2)):
            for k, z in enumerate((z1, z2)):
                radius = np.sqrt(x * x + y * y + z * z)
                term = _safe_xlog(x, y + radius) + _safe_xlog(y, x + radius)
                denominator = z * radius
                ratio = np.divide(
                    x * y,
                    denominator,
                    out=np.zeros_like(radius),
                    where=denominator != 0.0,
                )
                term -= z * np.arctan(ratio)
                total += (-1.0) ** (i + j + k) * term
    return G * np.asarray(density_kg_m3) * total


def nagy_prism_correction(
    cells: TerrainCells,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
) -> float:
    dh = cells.height_m - station_height_m
    x1 = cells.east_m - 0.5 * cells.width_east_m
    x2 = cells.east_m + 0.5 * cells.width_east_m
    y1 = cells.north_m - 0.5 * cells.width_north_m
    y2 = cells.north_m + 0.5 * cells.width_north_m
    z1 = np.minimum(0.0, dh)
    z2 = np.maximum(0.0, dh)
    contrast = np.sign(dh) * density_kg_m3
    acceleration = vertical_prism_acceleration(x1, x2, y1, y2, z1, z2, contrast)
    acceleration = np.where(np.isfinite(acceleration), acceleration, 0.0)
    return float(np.sum(acceleration, dtype=np.float64) * SI_TO_MGAL)
