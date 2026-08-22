"""Centroid point-mass approximation used only as a convergence diagnostic."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, G, SI_TO_MGAL
from ..io.terrain import TerrainCells


def point_mass_correction(
    cells: TerrainCells,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
) -> float:
    dh = cells.height_m - station_height_m
    half_height = 0.5 * dh
    radius_squared = cells.distance_m**2 + half_height**2
    valid = (radius_squared > 0.0) & (np.abs(dh) > 0.0)
    # Signed density contrast times signed centroid height makes both hills and
    # valleys a positive classical terrain correction.
    volume = cells.area_m2 * np.abs(dh)
    gz_up = G * density_kg_m3 * volume[valid] * np.abs(half_height[valid]) / (
        radius_squared[valid] ** 1.5
    )
    return float(np.sum(gz_up, dtype=np.float64) * SI_TO_MGAL)

