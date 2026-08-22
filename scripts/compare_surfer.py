#!/usr/bin/env python3
"""Compare original station heights with the legacy Surfer grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from terrain_correction_gravity.io.stations import load_stations
from terrain_correction_gravity.io.surfer import sample_xyz_grid


def metrics(error: np.ndarray) -> dict[str, float | int]:
    return {
        "n": int(error.size),
        "bias_m": float(error.mean()),
        "mae_m": float(np.abs(error).mean()),
        "rmse_m": float(np.sqrt(np.mean(error**2))),
        "median_absolute_error_m": float(np.median(np.abs(error))),
        "p95_absolute_error_m": float(np.quantile(np.abs(error), 0.95)),
        "maximum_absolute_error_m": float(np.abs(error).max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stations", type=Path, required=True)
    parser.add_argument("--out1", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stations = load_stations(args.stations)
    surfer = sample_xyz_grid(args.out1, stations.longitude_deg, stations.latitude_deg)
    error = surfer - stations.altitude_m
    table = pd.DataFrame(
        {
            "excel_row": stations.excel_row,
            "longitude_deg": stations.longitude_deg,
            "latitude_deg": stations.latitude_deg,
            "station_altitude_m": stations.altitude_m,
            "surfer_bilinear_altitude_m": surfer,
            "surfer_minus_station_m": error,
        }
    )
    table.to_csv(args.output_dir / "surfer_station_comparison.csv", index=False)
    positive = stations.altitude_m > 0.0
    report = {
        "definition": "error = bilinear out1.dat altitude - original station altitude",
        "all_stations": metrics(error),
        "positive_altitude_stations": metrics(error[positive]),
        "warning": (
            "out1.dat was interpolated from the same station heights; these are fidelity, "
            "not independent terrain-accuracy, metrics."
        ),
    }
    (args.output_dir / "surfer_comparison_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()

