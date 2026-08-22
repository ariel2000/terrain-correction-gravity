#!/usr/bin/env python3
"""Inventory DEM tiles and write portable integrity metadata."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fabdem-dir", type=Path, required=True)
    parser.add_argument("--copernicus-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for lat0 in range(1, 9):
        for lon0 in range(6, 13):
            tile = f"N{lat0:02d}E{lon0:03d}"
            products = {
                "FABDEM_v1.2": args.fabdem_dir / f"{tile}_FABDEM_V1-2.tif",
                "Copernicus_GLO30": args.copernicus_dir
                / f"Copernicus_DSM_COG_10_N{lat0:02d}_00_E{lon0:03d}_00_DEM.tif",
            }
            for product, path in products.items():
                present = path.exists()
                rows.append(
                    {
                        "product": product,
                        "tile": tile,
                        "filename": path.name,
                        "status": "present" if present else "absent_ocean_tile",
                        "size_bytes": path.stat().st_size if present else 0,
                        "sha256": sha256(path) if present else "",
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

