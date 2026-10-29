"""Original and station-centred spherical Vannes terrain-effect kernels."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, DEG_TO_RAD, EARTH_RADIUS, SI_TO_MGAL, G
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

    distance = np.sqrt(np.maximum(r * r + geoid_radius_m**2 - 2.0 * r * geoid_radius_m * x, 0.0))
    sigma = height * (1.0 + height / geoid_radius_m + height * height / (3.0 * geoid_radius_m**2))
    condensed = np.divide(
        -(geoid_radius_m**2) * sigma * (r - geoid_radius_m * x),
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
        density_kg_m3 * (source_kernel - reference_kernel) * cells.solid_angle_sr[nonself],
        dtype=np.float64,
    )
    return float(G * integral * SI_TO_MGAL)


def _select_cells(cells: TerrainCells, mask: np.ndarray) -> TerrainCells:
    """Return a cell subset without importing the higher-level pipeline module."""
    return TerrainCells(
        **{name: getattr(cells, name)[mask] for name in TerrainCells.__dataclass_fields__}
    )


def station_centered_spherical_kernel(
    computation_radius_m: float,
    cosine_psi: np.ndarray,
    terrain_radius_m: np.ndarray,
    quadrature_order: int = 8,
) -> np.ndarray:
    """Return the SCCV residual kernel per unit solid angle.

    The residual mass between the station sphere and the terrain is integrated
    in three dimensions and then condensed, with exactly conserved mass, onto
    the sphere passing through the station.  Positive and negative relief are
    represented by opposite density contrasts.  The returned kernel has units
    of metres and must be multiplied by ``G * density * solid_angle``.

    Cells at zero angular separation are intentionally excluded here.  The
    containing/near-field cells are evaluated analytically as Nagy prisms by
    :func:`vannes_adapted_correction`.
    """
    if quadrature_order < 2:
        raise ValueError("quadrature_order must be at least 2")
    a = float(computation_radius_m)
    x = np.clip(np.asarray(cosine_psi, dtype=np.float64), -1.0, 1.0)
    b = np.asarray(terrain_radius_m, dtype=np.float64)
    x, b = np.broadcast_arrays(x, b)
    dh = b - a
    sign = np.sign(dh)
    lower = np.minimum(a, b)
    upper = np.maximum(a, b)
    half = 0.5 * (upper - lower)
    midpoint = 0.5 * (upper + lower)

    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    radial = midpoint[:, None] + half[:, None] * nodes[None, :]
    distance_squared = np.maximum(a * a + radial * radial - 2.0 * a * radial * x[:, None], 0.0)
    distance = np.sqrt(distance_squared)
    numerator = radial * radial * (radial * x[:, None] - a)
    integrand = np.divide(
        numerator,
        distance**3,
        out=np.zeros_like(numerator),
        where=distance > 0.0,
    )
    topographic = sign * half * np.sum(integrand * weights[None, :], axis=1)

    # Signed volume per unit solid angle, condensed onto r=a.  This preserves
    # mass exactly, including the spherical r^2 dr volume factor.
    signed_volume = sign * (upper**3 - lower**3) / 3.0
    condensed_distance = a * np.sqrt(np.maximum(2.0 * (1.0 - x), 0.0))
    condensed = np.divide(
        signed_volume * a * (x - 1.0),
        condensed_distance**3,
        out=np.zeros_like(signed_volume),
        where=condensed_distance > 0.0,
    )
    return topographic - condensed


def vannes_adapted_correction(
    cells: TerrainCells,
    station_lon_deg: float,
    station_lat_deg: float,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
    geoid_radius_m: float = EARTH_RADIUS,
    near_field_radius_m: float = 5_000.0,
    quadrature_order: int = 8,
) -> float:
    """Compute the station-centred condensed-mass Vannes correction in mGal.

    Rectangular Nagy prisms are used only in the near field, where a spherical
    point-cell kernel is singular and the containing DEM cell must be treated
    analytically.  Outside that radius, the SCCV kernel uses spherical volume
    integration and station-level mass condensation.
    """
    if near_field_radius_m < 0.0:
        raise ValueError("near_field_radius_m must be non-negative")
    if len(cells) == 0:
        return 0.0

    near = cells.distance_m < near_field_radius_m
    total_mgal = 0.0
    if np.any(near):
        # Local import avoids a module-level dependency cycle and makes the
        # shared near-field reference explicit in the implementation.
        from .prism import nagy_prism_correction

        total_mgal += nagy_prism_correction(
            _select_cells(cells, near), station_height_m, density_kg_m3
        )

    far = ~near
    if not np.any(far):
        return float(total_mgal)
    far_cells = _select_cells(cells, far)
    source_lat = far_cells.latitude_deg * DEG_TO_RAD
    source_lon = far_cells.longitude_deg * DEG_TO_RAD
    station_lat = float(station_lat_deg) * DEG_TO_RAD
    station_lon = float(station_lon_deg) * DEG_TO_RAD
    cosine_psi = np.sin(source_lat) * np.sin(station_lat) + np.cos(source_lat) * np.cos(
        station_lat
    ) * np.cos(source_lon - station_lon)
    computation_radius = geoid_radius_m + float(station_height_m)
    terrain_radius = geoid_radius_m + far_cells.height_m
    valid = (np.abs(terrain_radius - computation_radius) > 0.0) & (cosine_psi < 1.0 - 1e-15)
    if not np.any(valid):
        return float(total_mgal)
    kernel = station_centered_spherical_kernel(
        computation_radius,
        cosine_psi[valid],
        terrain_radius[valid],
        quadrature_order=quadrature_order,
    )
    spherical = (
        G * density_kg_m3 * np.sum(kernel * far_cells.solid_angle_sr[valid], dtype=np.float64)
    )
    return float(total_mgal + spherical * SI_TO_MGAL)
