"""Readers for the legacy Surfer grids used only as a reference baseline."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def sample_xyz_grid(path: str | Path, longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("Expected an ASCII XYZ grid with at least three columns")
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    if x.size * y.size != data.shape[0]:
        raise ValueError("XYZ points do not form a complete regular grid")
    ordered = data[np.lexsort((data[:, 0], data[:, 1]))]
    z = ordered[:, 2].reshape(y.size, x.size)
    interpolator = RegularGridInterpolator((y, x), z, bounds_error=True)
    return interpolator(np.column_stack([latitude, longitude]))

