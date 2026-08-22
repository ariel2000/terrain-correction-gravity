import numpy as np

from terrain_correction_gravity.constants import G
from terrain_correction_gravity.methods.prism import vertical_prism_acceleration


def test_small_far_prism_converges_to_point_mass_above():
    density = 2670.0
    result = vertical_prism_acceleration(
        np.array([-0.5]),
        np.array([0.5]),
        np.array([-0.5]),
        np.array([0.5]),
        np.array([10.0]),
        np.array([11.0]),
        density,
    )[0]
    expected = G * density * 1.0 * 10.5 / 10.5**3
    assert np.isclose(result, expected, rtol=1e-4)


def test_prism_sign_reverses_below_observer():
    above = vertical_prism_acceleration(
        np.array([-0.5]), np.array([0.5]), np.array([-0.5]), np.array([0.5]),
        np.array([10.0]), np.array([11.0]), 2670.0,
    )[0]
    below = vertical_prism_acceleration(
        np.array([-0.5]), np.array([0.5]), np.array([-0.5]), np.array([0.5]),
        np.array([-11.0]), np.array([-10.0]), 2670.0,
    )[0]
    assert np.isclose(above, -below, rtol=1e-10)

