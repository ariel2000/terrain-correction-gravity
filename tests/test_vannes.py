import numpy as np

from terrain_correction_gravity.io.terrain import TerrainCells
from terrain_correction_gravity.methods.prism import nagy_prism_correction
from terrain_correction_gravity.methods.vannes import (
    dnewtc,
    station_centered_spherical_kernel,
    vannes_adapted_correction,
)


def test_dnewtc_vector_is_finite_away_from_self_point():
    x = np.cos(np.deg2rad(np.array([0.01, 0.1, 1.0])))
    result = dnewtc(6_371_500.0, x, 6_371_000.0, np.array([6_371_100.0] * 3), np.array([100.0] * 3))
    assert np.isfinite(result).all()


def _cells(height_m: float, distance_m: float = 10_000.0) -> TerrainCells:
    latitude = distance_m / 6_371_000.0 * 180.0 / np.pi
    width = 100.0
    solid_angle = width * width / 6_371_000.0**2
    values = np.array([0.0], dtype=float)
    return TerrainCells(
        longitude_deg=values.copy(),
        latitude_deg=np.array([latitude]),
        height_m=np.array([height_m]),
        east_m=values.copy(),
        north_m=np.array([distance_m]),
        width_east_m=np.array([width]),
        width_north_m=np.array([width]),
        area_m2=np.array([width * width]),
        distance_m=np.array([distance_m]),
        azimuth_rad=values.copy(),
        solid_angle_sr=np.array([solid_angle]),
    )


def test_station_centered_kernel_is_zero_for_flat_terrain():
    result = station_centered_spherical_kernel(
        6_371_500.0, np.array([0.99]), np.array([6_371_500.0])
    )
    assert result[0] == 0.0


def test_adapted_correction_is_positive_for_hill_and_valley():
    hill = vannes_adapted_correction(_cells(700.0), 0.0, 0.0, 500.0, near_field_radius_m=0.0)
    valley = vannes_adapted_correction(_cells(300.0), 0.0, 0.0, 500.0, near_field_radius_m=0.0)
    assert hill > 0.0
    assert valley > 0.0


def test_adapted_correction_is_linear_in_density():
    cells = _cells(800.0)
    first = vannes_adapted_correction(cells, 0.0, 0.0, 500.0, density_kg_m3=1000.0)
    second = vannes_adapted_correction(cells, 0.0, 0.0, 500.0, density_kg_m3=2000.0)
    assert np.isclose(second, 2.0 * first, rtol=1e-12, atol=1e-15)


def test_quadrature_orders_converge():
    cells = _cells(900.0, distance_m=20_000.0)
    order8 = vannes_adapted_correction(
        cells, 0.0, 0.0, 500.0, near_field_radius_m=0.0, quadrature_order=8
    )
    order16 = vannes_adapted_correction(
        cells, 0.0, 0.0, 500.0, near_field_radius_m=0.0, quadrature_order=16
    )
    assert np.isclose(order8, order16, rtol=1e-8, atol=1e-12)


def test_near_field_branch_equals_nagy():
    cells = _cells(900.0, distance_m=1_000.0)
    adapted = vannes_adapted_correction(cells, 0.0, 0.0, 500.0, near_field_radius_m=5_000.0)
    reference = nagy_prism_correction(cells, 500.0)
    assert adapted == reference
