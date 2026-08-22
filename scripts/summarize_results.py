#!/usr/bin/env python3
"""Create compact tables, figures and a Markdown scientific summary."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHODS = {
    "hammer_mgal": "Hammer-style",
    "nagy_prism_mgal": "Nagy prism",
    "point_mass_mgal": "Point mass",
    "vannes_direct_mgal": "Vannes direct",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    figures.mkdir(exist_ok=True)
    data = pd.read_csv(args.input_csv)
    if not np.isfinite(data[list(METHODS)].to_numpy()).all():
        raise ValueError("Non-finite terrain result detected")
    expected_rows = data["excel_row"].nunique() * data["radius_km"].nunique()
    if len(data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} unique station-radius rows, found {len(data)}")
    shutil.copyfile(args.input_csv, args.output_dir / "terrain_corrections.csv")

    summaries = []
    for radius, group in data.groupby("radius_km"):
        for column, label in METHODS.items():
            values = group[column]
            summaries.append(
                {
                    "radius_km": radius,
                    "method": label,
                    "n": len(values),
                    "mean_mgal": values.mean(),
                    "std_mgal": values.std(ddof=1),
                    "minimum_mgal": values.min(),
                    "median_mgal": values.median(),
                    "p95_mgal": values.quantile(0.95),
                    "maximum_mgal": values.max(),
                    "negative_count": int((values < 0).sum()),
                }
            )
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output_dir / "summary_statistics.csv", index=False)

    radii = sorted(data["radius_km"].unique())
    wide = data.pivot(index="excel_row", columns="radius_km", values=list(METHODS))
    convergence_rows = []
    for column, label in METHODS.items():
        for inner, outer in zip(radii[:-1], radii[1:]):
            increment = wide[column][outer] - wide[column][inner]
            convergence_rows.append(
                {
                    "method": label,
                    "inner_radius_km": inner,
                    "outer_radius_km": outer,
                    "mean_increment_mgal": increment.mean(),
                    "median_absolute_increment_mgal": increment.abs().median(),
                    "p95_absolute_increment_mgal": increment.abs().quantile(0.95),
                    "maximum_absolute_increment_mgal": increment.abs().max(),
                    "fraction_below_0_01_mgal": (increment.abs() < 0.01).mean(),
                }
            )
    convergence = pd.DataFrame(convergence_rows)
    convergence.to_csv(args.output_dir / "convergence_statistics.csv", index=False)

    final = data[data["radius_km"] == max(radii)].copy()
    final.to_csv(args.output_dir / "station_results_200km.csv", index=False)
    qc = (
        final.groupby("station_qc_status", dropna=False)
        .agg(stations=("excel_row", "count"), maximum_prism_mgal=("nagy_prism_mgal", "max"))
        .reset_index()
    )
    qc.to_csv(args.output_dir / "quality_control_summary.csv", index=False)

    plt.figure(figsize=(9, 5.5))
    plt.boxplot(
        [final[column] for column in METHODS],
        tick_labels=list(METHODS.values()),
        showfliers=False,
    )
    plt.ylabel("Terrain effect / correction (mGal)")
    plt.title("Method distributions at 200 km (outliers hidden in display only)")
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(figures / "method_distributions_200km.png", dpi=180)
    plt.close()

    plt.figure(figsize=(9, 5.5))
    for label, group in convergence.groupby("method"):
        x = [f"{int(a)}-{int(b)}" for a, b in zip(group.inner_radius_km, group.outer_radius_km)]
        plt.plot(x, group.p95_absolute_increment_mgal, marker="o", label=label)
    plt.axhline(0.01, color="black", linestyle="--", linewidth=1, label="0.01 mGal criterion")
    plt.yscale("log")
    plt.ylabel("95th percentile of |increment| (mGal)")
    plt.xlabel("Added radial interval (km)")
    plt.title("Radial convergence remains incomplete at 200 km")
    plt.grid(alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures / "radial_convergence.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7, 7))
    scatter = plt.scatter(
        final.longitude_deg,
        final.latitude_deg,
        c=final.nagy_prism_mgal,
        s=13,
        cmap="viridis",
        vmin=0,
        vmax=final.nagy_prism_mgal.quantile(0.98),
    )
    plt.colorbar(scatter, label="Nagy-prism correction (mGal)")
    plt.xlabel("Longitude (deg E)")
    plt.ylabel("Latitude (deg N)")
    plt.title("Terrain correction at 200 km")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(figures / "map_nagy_prism_200km.png", dpi=180)
    plt.close()

    final_summary = summary[summary.radius_km == max(radii)].set_index("method")
    last_increment = convergence[
        (convergence.inner_radius_km == radii[-2])
        & (convergence.outer_radius_km == radii[-1])
    ].set_index("method")
    markdown = f"""# Results

The complete run contains {len(data):,} station-radius rows ({final.excel_row.nunique():,}
stations and {len(radii)} radii). All reported method values are finite.

## Distribution at {max(radii):g} km

| Method | Mean | Median | P95 | Maximum | Negative stations |
|---|---:|---:|---:|---:|---:|
"""
    for label in METHODS.values():
        row = final_summary.loc[label]
        markdown += (
            f"| {label} | {row.mean_mgal:.3f} | {row.median_mgal:.3f} | "
            f"{row.p95_mgal:.3f} | {row.maximum_mgal:.3f} | {int(row.negative_count)} |\n"
        )
    markdown += "\n## Convergence from 100 to 200 km\n\n"
    markdown += "| Method | Mean increment | P95 absolute increment | Fraction below 0.01 mGal |\n"
    markdown += "|---|---:|---:|---:|\n"
    for label in METHODS.values():
        row = last_increment.loc[label]
        markdown += (
            f"| {label} | {row.mean_increment_mgal:.3f} | "
            f"{row.p95_absolute_increment_mgal:.3f} | "
            f"{100*row.fraction_below_0_01_mgal:.1f}% |\n"
        )
    markdown += """

The 0.01 mGal convergence target is not reached by most stations at 200 km.
The 200 km estimates are therefore the largest-domain results supported by the
downloaded DEM, not mathematically converged infinite-domain corrections.

Nagy-prism is the reference conventional correction because it uses the exact
rectangular-cell solution. Hammer-style is lower on average because angular and
radial compartment averaging removes within-compartment relief variance. The
point-mass estimate agrees closely with Nagy in the far field but differs in the
near field. Vannes remains signed and can be negative because it is a spherical
residual direct effect, not the always-positive conventional correction.

Stations excluded from DEM conditioning remain in the gravity results. Their
terrain corrections, especially the largest values, require independent checking
of station elevation and vertical datum before geological interpretation.
"""
    (args.output_dir / "RESULTS.md").write_text(markdown, encoding="utf-8")


if __name__ == "__main__":
    main()

