"""Command-line interface for audit and calculation."""

from __future__ import annotations

import argparse
import csv
import json
import tomllib
from pathlib import Path

from .io.stations import audit_stations, load_stations
from .io.terrain import TiledTerrain
from .pipeline import ResolutionTier, calculate_station_radii


def _config(path: Path) -> dict:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def _resolved(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def command_audit(args: argparse.Namespace) -> None:
    stations = load_stations(args.stations)
    print(json.dumps(audit_stations(stations), indent=2, ensure_ascii=False))


def command_run(args: argparse.Namespace) -> None:
    cfg_path = Path(args.config).resolve()
    cfg = _config(cfg_path)
    root = cfg_path.parent.parent
    paths = cfg["paths"]
    stations = load_stations(_resolved(root, paths["stations"]))
    excluded_by_row: dict[int, list[str]] = {}
    exclusions_path = _resolved(root, paths["excluded_controls"])
    if exclusions_path.exists():
        exclusions = json.loads(exclusions_path.read_text(encoding="utf-8"))
        for control in exclusions:
            for row in control["source_rows"]:
                excluded_by_row[int(row)] = list(control["exclusion_reasons"])
    terrain = TiledTerrain(
        _resolved(root, paths["hybrid_dir"]),
        _resolved(root, paths["fabdem_dir"]),
    )
    tiers = [ResolutionTier(**item) for item in cfg["resolution_tiers"]]
    radii = args.radii or cfg["domain"]["radii_km"]
    output_dir = _resolved(root, paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "terrain_corrections.csv"
    fieldnames = [
        "excel_row",
        "longitude_deg",
        "latitude_deg",
        "station_altitude_m",
        "free_air_anomaly_mgal",
        "station_qc_status",
        "station_qc_reasons",
        "radius_km",
        "cell_count",
        "hammer_mgal",
        "nagy_prism_mgal",
        "point_mass_mgal",
        "vannes_direct_mgal",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(len(stations)):
            station_results = calculate_station_radii(
                terrain,
                float(stations.longitude_deg[index]),
                float(stations.latitude_deg[index]),
                float(stations.altitude_m[index]),
                [float(radius) for radius in radii],
                tiers,
                float(cfg["physics"]["density_kg_m3"]),
                int(cfg["hammer"]["sectors"]),
                int(cfg["hammer"]["rings_per_decade"]),
            )
            for radius, result in zip(radii, station_results, strict=True):
                writer.writerow(
                    {
                        "excel_row": int(stations.excel_row[index]),
                        "longitude_deg": stations.longitude_deg[index],
                        "latitude_deg": stations.latitude_deg[index],
                        "station_altitude_m": stations.altitude_m[index],
                        "free_air_anomaly_mgal": stations.free_air_anomaly_mgal[index],
                        "station_qc_status": (
                            "excluded_from_dem_conditioning"
                            if int(stations.excel_row[index]) in excluded_by_row
                            else "accepted_for_dem_conditioning"
                        ),
                        "station_qc_reasons": ";".join(
                            excluded_by_row.get(int(stations.excel_row[index]), [])
                        ),
                        "radius_km": radius,
                        **result,
                    }
                )
                stream.flush()
            if (index + 1) % 25 == 0 or index + 1 == len(stations):
                print(f"Completed {index + 1}/{len(stations)} stations", flush=True)
    print(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(prog="terrain-correction")
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit-stations")
    audit.add_argument("stations", type=Path)
    audit.set_defaults(func=command_audit)
    run = subparsers.add_parser("run")
    run.add_argument("--config", default="config/default.toml")
    run.add_argument("--radii", nargs="+", type=float)
    run.set_defaults(func=command_run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
