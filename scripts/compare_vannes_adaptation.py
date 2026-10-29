"""Controlled SCCV comparison on the information-equivalent historical grid."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from terrain_correction_gravity.io.stations import load_stations
from terrain_correction_gravity.io.terrain import TerrainCells
from terrain_correction_gravity.methods.hammer import (
    hammer_sector_correction,
    logarithmic_edges,
)
from terrain_correction_gravity.methods.point_mass import point_mass_correction
from terrain_correction_gravity.methods.prism import nagy_prism_correction
from terrain_correction_gravity.methods.vannes import (
    vannes_adapted_correction,
    vannes_direct_effect,
)

ROOT = Path(__file__).resolve().parents[1]
EARTH_RADIUS_M = 6_371_000.0
DEG_TO_RAD = np.pi / 180.0

_GRID: dict[str, np.ndarray] = {}


def load_information_equivalent_grid(path: Path) -> dict[str, np.ndarray]:
    """Reduce the 598 x 598 Surfer interpolation to its 200 x 200 support."""
    data = np.loadtxt(path)
    longitude = np.unique(data[:, 0])
    latitude = np.unique(data[:, 1])
    if longitude.size != 598 or latitude.size != 598:
        raise ValueError("Expected the historical 598 x 598 grid")
    height = data[:, 2].reshape(latitude.size, longitude.size)
    longitude = longitude[::3]
    latitude = latitude[::3]
    height = height[::3, ::3]
    if height.shape != (200, 200):
        raise ValueError("The information-equivalent grid must be 200 x 200")
    lon_grid, lat_grid = np.meshgrid(longitude, latitude)
    step_lon_rad = float(np.median(np.diff(longitude))) * DEG_TO_RAD
    step_lat_rad = float(np.median(np.diff(latitude))) * DEG_TO_RAD
    lat_rad = lat_grid * DEG_TO_RAD
    width_north = np.full_like(lat_grid, EARTH_RADIUS_M * step_lat_rad)
    width_east = EARTH_RADIUS_M * np.cos(lat_rad) * step_lon_rad
    solid_angle = step_lon_rad * (
        np.sin(lat_rad + 0.5 * step_lat_rad) - np.sin(lat_rad - 0.5 * step_lat_rad)
    )
    return {
        "longitude_deg": lon_grid.ravel().astype(np.float64),
        "latitude_deg": lat_grid.ravel().astype(np.float64),
        "height_m": height.ravel().astype(np.float64),
        "width_east_m": width_east.ravel().astype(np.float64),
        "width_north_m": width_north.ravel().astype(np.float64),
        "solid_angle_sr": solid_angle.ravel().astype(np.float64),
        "area_m2": (EARTH_RADIUS_M**2 * solid_angle).ravel().astype(np.float64),
    }


def _initialise_worker(grid: dict[str, np.ndarray]) -> None:
    global _GRID
    _GRID = grid


def station_cells(longitude_deg: float, latitude_deg: float) -> TerrainCells:
    dlat = (_GRID["latitude_deg"] - latitude_deg) * DEG_TO_RAD
    dlon = (_GRID["longitude_deg"] - longitude_deg) * DEG_TO_RAD
    mean_lat = 0.5 * (_GRID["latitude_deg"] + latitude_deg) * DEG_TO_RAD
    north = EARTH_RADIUS_M * dlat
    east = EARTH_RADIUS_M * np.cos(mean_lat) * dlon
    distance = np.hypot(east, north)
    return TerrainCells(
        longitude_deg=_GRID["longitude_deg"],
        latitude_deg=_GRID["latitude_deg"],
        height_m=_GRID["height_m"],
        east_m=east,
        north_m=north,
        width_east_m=_GRID["width_east_m"],
        width_north_m=_GRID["width_north_m"],
        area_m2=_GRID["area_m2"],
        distance_m=distance,
        azimuth_rad=np.mod(np.arctan2(east, north), 2.0 * np.pi),
        solid_angle_sr=_GRID["solid_angle_sr"],
    )


def calculate_one(item: tuple[int, int, float, float, float, float]) -> dict[str, object]:
    index, excel_row, longitude, latitude, altitude, free_air = item
    cells = station_cells(longitude, latitude)
    maximum_radius = max(float(np.max(cells.distance_m)), 1.0)
    hammer_edges = logarithmic_edges(maximum_radius * 1.000001)
    nagy = nagy_prism_correction(cells, altitude)
    adapted = vannes_adapted_correction(
        cells,
        longitude,
        latitude,
        altitude,
        near_field_radius_m=5_000.0,
        quadrature_order=8,
    )
    original = vannes_direct_effect(cells, longitude, latitude, altitude)
    return {
        "station_id": f"CMG-S{index + 1:04d}",
        "excel_row": excel_row,
        "longitude_deg": longitude,
        "latitude_deg": latitude,
        "station_altitude_m": altitude,
        "free_air_anomaly_mgal": free_air,
        "cell_count": len(cells),
        "domain": "historical_grid_200x200_full_rectangle",
        "near_field_radius_km": 5.0,
        "hammer_mgal": hammer_sector_correction(cells, altitude, hammer_edges, 12),
        "nagy_prism_mgal": nagy,
        "point_mass_mgal": point_mass_correction(cells, altitude),
        "vannes_direct_mgal": original,
        "vannes_adapted_mgal": adapted,
        "sccv_minus_nagy_mgal": adapted - nagy,
        "abs_sccv_minus_nagy_mgal": abs(adapted - nagy),
    }


def write_summary(rows: list[dict[str, object]], path: Path) -> None:
    nagy = np.array([float(row["nagy_prism_mgal"]) for row in rows])
    methods = {
        "Vannes original": np.array([float(row["vannes_direct_mgal"]) for row in rows]),
        "SCCV": np.array([float(row["vannes_adapted_mgal"]) for row in rows]),
        "Hammer-style": np.array([float(row["hammer_mgal"]) for row in rows]),
        "Point mass": np.array([float(row["point_mass_mgal"]) for row in rows]),
    }
    with path.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = [
            "method",
            "n",
            "mean_bias_mgal",
            "mae_mgal",
            "rmse_mgal",
            "p95_abs_error_mgal",
            "maximum_abs_error_mgal",
            "pearson_r",
            "negative_count",
            "within_0_01_mgal_percent",
            "within_0_05_mgal_percent",
            "within_0_10_mgal_percent",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for name, values in methods.items():
            error = values - nagy
            absolute = np.abs(error)
            writer.writerow(
                {
                    "method": name,
                    "n": values.size,
                    "mean_bias_mgal": float(np.mean(error)),
                    "mae_mgal": float(np.mean(absolute)),
                    "rmse_mgal": float(np.sqrt(np.mean(error**2))),
                    "p95_abs_error_mgal": float(np.percentile(absolute, 95)),
                    "maximum_abs_error_mgal": float(np.max(absolute)),
                    "pearson_r": float(np.corrcoef(values, nagy)[0, 1]),
                    "negative_count": int(np.count_nonzero(values < 0.0)),
                    "within_0_01_mgal_percent": float(np.mean(absolute <= 0.01) * 100.0),
                    "within_0_05_mgal_percent": float(np.mean(absolute <= 0.05) * 100.0),
                    "within_0_10_mgal_percent": float(np.mean(absolute <= 0.10) * 100.0),
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-stations", type=int)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "outputs" / "vannes_adapted_reference_grid"
    )
    args = parser.parse_args()
    grid = load_information_equivalent_grid(ROOT / "data" / "reference" / "out1_surfer.dat")
    stations = load_stations(ROOT / "data" / "raw" / "original_stations.xlsx")
    count = len(stations) if args.max_stations is None else min(args.max_stations, len(stations))
    items = [
        (
            index,
            int(stations.excel_row[index]),
            float(stations.longitude_deg[index]),
            float(stations.latitude_deg[index]),
            float(stations.altitude_m[index]),
            float(stations.free_air_anomaly_mgal[index]),
        )
        for index in range(count)
    ]
    if args.workers == 1:
        _initialise_worker(grid)
        rows = [calculate_one(item) for item in items]
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers, initializer=_initialise_worker, initargs=(grid,)
        ) as pool:
            rows = list(pool.map(calculate_one, items, chunksize=4))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    station_path = args.output_dir / "station_comparison.csv"
    with station_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_summary(rows, args.output_dir / "pairwise_summary.csv")
    print(station_path)


if __name__ == "__main__":
    main()
