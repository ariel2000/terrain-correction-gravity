"""Spherical direct topographic effect transcribed from Vannes (2011), DNEWTC."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, DEG_TO_RAD, EARTH_RADIUS, G, SI_TO_MGAL
from ..io.terrain import TerrainCells


def dnewtc(
    computation_radius_m: float,
    cosine_psi: np.ndarray,
    lower_radius_m: float,
    upper_radius_m: np.ndarray,
    height_m: np.ndarray,
    geoid_radius_m: float = EARTH_RADIUS,
) -> np.ndarray:
    """First derivative of the residual Newton kernel from Vannes Appendix A.7.2."""
    r = float(computation_radius_m)
    x = np.asarray(cosine_psi, dtype=np.float64)
    lower = float(lower_radius_m)
    upper = np.asarray(upper_radius_m, dtype=np.float64)
    height = np.asarray(height_m, dtype=np.float64)
    t1 = lower / r
    t2 = upper / r
    d1 = np.sqrt(np.maximum(1.0 - 2.0 * t1 * x + t1 * t1, 0.0))
    d2 = np.sqrt(np.maximum(1.0 - 2.0 * t2 * x + t2 * t2, 0.0))
    numerator = t2 - x + d2
    denominator = t1 - x + d1
    ratio = np.divide(
        numerator,
        denominator,
        out=np.ones_like(numerator),
        where=np.abs(denominator) > np.finfo(float).eps,
    )
    log_term = (3.0 * x * x - 1.0) * np.log(np.maximum(np.abs(ratio), np.finfo(float).tiny))
    p1 = (x * t1 * t1 + 3.0 * x + (1.0 - 6.0 * x * x) * t1) / d1
    p2 = (x * t2 * t2 + 3.0 * x + (1.0 - 6.0 * x * x) * t2) / d2
    topographic = r * (p2 - p1 + log_term)

    distance = np.sqrt(
        np.maximum(r * r + geoid_radius_m**2 - 2.0 * r * geoid_radius_m * x, 0.0)
    )
    sigma = height * (
        1.0 + height / geoid_radius_m + height * height / (3.0 * geoid_radius_m**2)
    )
    condensed = np.divide(
        -geoid_radius_m**2 * sigma * (r - geoid_radius_m * x),
        distance**3,
        out=np.zeros_like(distance),
        where=distance > 0.0,
    )
    return topographic - condensed


def vannes_direct_effect(
    cells: TerrainCells,
    station_lon_deg: float,
    station_lat_deg: float,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
    geoid_radius_m: float = EARTH_RADIUS,
) -> float:
    """Compute Vannes' signed direct residual effect on surface gravity in mGal."""
    source_height = np.maximum(cells.height_m, 0.0)
    station_height = max(float(station_height_m), 0.0)
    source_lat = cells.latitude_deg * DEG_TO_RAD
    source_lon = cells.longitude_deg * DEG_TO_RAD
    station_lat = station_lat_deg * DEG_TO_RAD
    station_lon = station_lon_deg * DEG_TO_RAD
    cosine_psi = np.sin(source_lat) * np.sin(station_lat) + np.cos(source_lat) * np.cos(
        station_lat
    ) * np.cos(source_lon - station_lon)
    cosine_psi = np.clip(cosine_psi, -1.0, 1.0)
    nonself = cosine_psi < 1.0 - 1e-15
    if not np.any(nonself):
        return 0.0
    x = cosine_psi[nonself]
    r = geoid_radius_m + station_height
    source_kernel = dnewtc(
        r,
        x,
        geoid_radius_m,
        geoid_radius_m + source_height[nonself],
        source_height[nonself],
        geoid_radius_m,
    )
    reference_height = np.full_like(x, station_height)
    reference_kernel = dnewtc(
        r,
        x,
        geoid_radius_m,
        geoid_radius_m + reference_height,
        reference_height,
        geoid_radius_m,
    )
    integral = np.sum(
        density_kg_m3
        * (source_kernel - reference_kernel)
        * cells.solid_angle_sr[nonself],
        dtype=np.float64,
    )
    return float(G * integral * SI_TO_MGAL)

