import numpy as np

from terrain_correction_gravity.methods.vannes import dnewtc


def test_dnewtc_vector_is_finite_away_from_self_point():
    x = np.cos(np.deg2rad(np.array([0.01, 0.1, 1.0])))
    result = dnewtc(6_371_500.0, x, 6_371_000.0, np.array([6_371_100.0] * 3), np.array([100.0] * 3))
    assert np.isfinite(result).all()

