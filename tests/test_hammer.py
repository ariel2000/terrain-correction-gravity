import numpy as np

from terrain_correction_gravity.io.terrain import TerrainCells
from terrain_correction_gravity.methods.hammer import hammer_sector_correction, logarithmic_edges


def cells_with_height(height):
    one = np.array([1.0])
    return TerrainCells(
        longitude_deg=np.array([0.0]),
        latitude_deg=np.array([0.0]),
        height_m=np.array([height]),
        east_m=one,
        north_m=one,
        width_east_m=one,
        width_north_m=one,
        area_m2=one,
        distance_m=np.array([50.0]),
        azimuth_rad=np.array([0.0]),
        solid_angle_sr=np.array([1e-12]),
    )


def test_flat_terrain_has_zero_hammer_correction():
    result = hammer_sector_correction(cells_with_height(100.0), 100.0, np.array([0.0, 100.0]))
    assert abs(result) < 1e-15


def test_hill_and_valley_have_same_hammer_magnitude():
    hill = hammer_sector_correction(cells_with_height(120.0), 100.0, np.array([0.0, 100.0]))
    valley = hammer_sector_correction(cells_with_height(80.0), 100.0, np.array([0.0, 100.0]))
    assert np.isclose(hill, valley)
    assert hill > 0.0


def test_logarithmic_ring_boundaries_are_nested_across_radii():
    short = logarithmic_edges(20_000.0)
    long = logarithmic_edges(200_000.0)
    # Every complete boundary below 20 km is unchanged when the domain grows.
    assert np.allclose(short[:-1], long[: short.size - 1])
