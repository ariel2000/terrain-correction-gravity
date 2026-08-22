"""Hammer-style annular-sector approximation of the classical terrain correction."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, SI_TO_MGAL, G
from ..io.terrain import TerrainCells


def hammer_sector_correction(
    cells: TerrainCells,
    station_height_m: float,
    radial_edges_m: np.ndarray,
    sectors: int = 12,
    density_kg_m3: float = DEFAULT_DENSITY,
) -> float:
    """Aggregate DEM cells into compartments and apply Hammer's sector formula."""
    edges = np.asarray(radial_edges_m, dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("radial_edges_m must be a strictly increasing vector")
    if sectors < 1:
        raise ValueError("sectors must be positive")
    ring = np.searchsorted(edges, cells.distance_m, side="right") - 1
    sector_width = 2.0 * np.pi / sectors
    sector = np.floor(cells.azimuth_rad / sector_width).astype(int)
    valid = (ring >= 0) & (ring < edges.size - 1) & (sector >= 0) & (sector < sectors)
    if not np.any(valid):
        return 0.0
    group = ring[valid] * sectors + sector[valid]
    count = (edges.size - 1) * sectors
    weights = cells.area_m2[valid]
    weighted_height = np.bincount(
        group, weights=weights * cells.height_m[valid], minlength=count
    )
    area = np.bincount(group, weights=weights, minlength=count)
    occupied = area > 0.0
    mean_height = np.zeros(count, dtype=np.float64)
    mean_height[occupied] = weighted_height[occupied] / area[occupied]
    ring_index = np.arange(count) // sectors
    r1 = edges[ring_index]
    r2 = edges[ring_index + 1]
    dh = mean_height - station_height_m
    compartment = G * density_kg_m3 * sector_width * (
        (r2 - r1) + np.sqrt(r1**2 + dh**2) - np.sqrt(r2**2 + dh**2)
    )
    return float(np.sum(compartment[occupied], dtype=np.float64) * SI_TO_MGAL)


def logarithmic_edges(
    maximum_radius_m: float,
    minimum_radius_m: float = 15.0,
    rings_per_decade: int = 8,
) -> np.ndarray:
    if maximum_radius_m <= minimum_radius_m:
        return np.array([0.0, maximum_radius_m], dtype=np.float64)
    ratio = 10.0 ** (1.0 / rings_per_decade)
    positive = [minimum_radius_m]
    while positive[-1] * ratio < maximum_radius_m:
        positive.append(positive[-1] * ratio)
    positive.append(maximum_radius_m)
    positive = np.asarray(positive, dtype=np.float64)
    return np.concatenate([[0.0], positive])
