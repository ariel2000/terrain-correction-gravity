#!/usr/bin/env python3
"""Generate the scientific maps and comparison charts published in the repository.

The committed PNG files are derived products. Recreate them by pointing this
script at the twelve 1-degree hybrid tiles and the corresponding FABDEM tiles.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LightSource, TwoSlopeNorm
from PIL import Image


TILE_RE = re.compile(r"N(?P<lat>\d{2})E(?P<lon>\d{3})")
EXTENT = (8.0, 11.0, 3.0, 7.0)
STUDY_EXTENT = (8.502, 10.9933, 3.5006, 6.5)


def tile_key(path: Path) -> tuple[int, int]:
    match = TILE_RE.search(path.name)
    if not match:
        raise ValueError(f"Cannot infer tile coordinates from {path.name}")
    return int(match.group("lat")), int(match.group("lon"))


def read_tile(path: Path, pixels: int) -> np.ndarray:
    with Image.open(path) as image:
        image = image.convert("F").resize((pixels, pixels), Image.Resampling.BOX)
        return np.asarray(image, dtype=np.float32)


def mosaic(directory: Path, pattern: str, pixels: int = 360) -> np.ndarray:
    result = np.full((4 * pixels, 3 * pixels), np.nan, dtype=np.float32)
    paths = sorted(directory.glob(pattern))
    if len(paths) < 12:
        raise FileNotFoundError(f"Expected 12 tiles in {directory}, found {len(paths)}")
    for path in paths:
        lat, lon = tile_key(path)
        if not (3 <= lat <= 6 and 8 <= lon <= 10):
            continue
        row = 6 - lat
        col = lon - 8
        result[
            row * pixels : (row + 1) * pixels,
            col * pixels : (col + 1) * pixels,
        ] = read_tile(path, pixels)
    return result


def relief_rgb(elevation: np.ndarray, vmax: float = 2400.0) -> np.ndarray:
    terrain = np.ma.masked_invalid(np.clip(elevation, 0, vmax))
    shade = LightSource(azdeg=315, altdeg=40).shade(
        terrain,
        cmap=plt.get_cmap("terrain"),
        vert_exag=1.2,
        blend_mode="soft",
        vmin=0,
        vmax=vmax,
    )
    shade[..., 3] = np.where(np.isfinite(elevation), 1.0, 0.0)
    return shade


def load_boundary(path: Path) -> np.ndarray:
    feature = json.loads(path.read_text(encoding="utf-8"))["features"][0]
    return np.asarray(feature["geometry"]["coordinates"][0], dtype=float)


def decorate_map(ax, title: str) -> None:
    ax.set_title(title, weight="bold")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="0.3", alpha=0.18, linewidth=0.6)


def study_area_figure(
    hybrid: np.ndarray,
    stations: pd.DataFrame,
    boundary: np.ndarray,
    output: Path,
) -> None:
    fig = plt.figure(figsize=(13, 8.4), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=(1.15, 3.2))
    inset = fig.add_subplot(grid[0, 0])
    inset.fill(boundary[:, 0], boundary[:, 1], color="#e6e0d4", edgecolor="#3f3f3f", linewidth=1.1)
    west, east, south, north = STUDY_EXTENT
    inset.add_patch(
        plt.Rectangle(
            (west, south),
            east - west,
            north - south,
            fill=False,
            edgecolor="#b2182b",
            linewidth=2.0,
        )
    )
    inset.text(11.7, 7.5, "CAMEROUN", ha="center", va="center", fontsize=11, weight="bold")
    inset.set_xlim(7.7, 16.5)
    inset.set_ylim(1.2, 13.2)
    decorate_map(inset, "Localisation nationale")

    ax = fig.add_subplot(grid[0, 1])
    ax.imshow(relief_rgb(hybrid), extent=EXTENT, origin="upper", interpolation="bilinear")
    points = ax.scatter(
        stations.longitude_deg,
        stations.latitude_deg,
        c=stations.station_altitude_m,
        s=7,
        cmap="magma",
        vmin=0,
        vmax=np.nanpercentile(stations.station_altitude_m, 98),
        edgecolors="none",
        alpha=0.82,
    )
    ax.plot(boundary[:, 0], boundary[:, 1], color="white", linewidth=1.0, alpha=0.9)
    ax.set_xlim(8.35, 11.15)
    ax.set_ylim(3.35, 6.65)
    decorate_map(ax, "MNT hybride à 30 m et 1 991 stations gravimétriques")
    cbar = fig.colorbar(points, ax=ax, shrink=0.78, pad=0.02)
    cbar.set_label("Altitude de la station (m)")
    fig.suptitle("Zone d’étude : Sud-Ouest et Ouest du Cameroun", fontsize=16, weight="bold")
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def dem_comparison_figure(
    fabdem: np.ndarray,
    hybrid: np.ndarray,
    surfer_path: Path,
    stations: pd.DataFrame,
    output: Path,
) -> None:
    raw = np.loadtxt(surfer_path)
    nx = np.unique(raw[:, 0]).size
    ny = np.unique(raw[:, 1]).size
    surfer = raw[:, 2].reshape(ny, nx)
    surfer = np.flipud(surfer)
    px = hybrid.shape[1] / (EXTENT[1] - EXTENT[0])
    py = hybrid.shape[0] / (EXTENT[3] - EXTENT[2])
    left = round((STUDY_EXTENT[0] - EXTENT[0]) * px)
    right = round((STUDY_EXTENT[1] - EXTENT[0]) * px)
    top = round((EXTENT[3] - STUDY_EXTENT[3]) * py)
    bottom = round((EXTENT[3] - STUDY_EXTENT[2]) * py)
    crop = hybrid[top:bottom, left:right]
    with Image.fromarray(crop.astype(np.float32), mode="F") as image:
        resized = image.resize((nx, ny), Image.Resampling.BILINEAR)
        hybrid_on_surfer = np.asarray(resized, dtype=float)
    difference = hybrid_on_surfer - surfer

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16, 5.8),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
    datasets = [surfer, hybrid_on_surfer, difference]
    titles = ["Ancien MNT Surfer", "MNT hybride FABDEM–stations", "Hybride − Surfer"]
    cmaps = ["terrain", "terrain", "RdBu_r"]
    norms = [None, None, TwoSlopeNorm(vcenter=0, vmin=-300, vmax=300)]
    for ax, values, title, cmap, norm in zip(axes, datasets, titles, cmaps, norms):
        image = ax.imshow(
            values,
            extent=STUDY_EXTENT,
            origin="upper",
            cmap=cmap,
            vmin=None if norm else 0,
            vmax=None if norm else 2300,
            norm=norm,
            interpolation="bilinear",
        )
        ax.scatter(stations.longitude_deg, stations.latitude_deg, s=1.2, color="black", alpha=0.25)
        decorate_map(ax, title)
        cbar = fig.colorbar(image, ax=ax, shrink=0.75, pad=0.015)
        cbar.set_label("Écart d’altitude (m)" if norm else "Altitude (m)")
    fig.suptitle(
        "Comparaison spatiale des modèles numériques de terrain",
        fontsize=15,
        weight="bold",
    )
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def validation_figure(report: dict, output: Path) -> None:
    cv = report["spatial_cross_validation"]
    base = cv["fabdem_baseline"]
    hybrid = cv["hybrid_cross_validated"]
    labels = ["MAE", "RMSE", "P95 |erreur|"]
    base_values = [base["mae_m"], base["rmse_m"], base["p95_absolute_error_m"]]
    hybrid_values = [hybrid["mae_m"], hybrid["rmse_m"], hybrid["p95_absolute_error_m"]]
    y = np.arange(len(labels))
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.5, 5.2),
        gridspec_kw={"width_ratios": (2.5, 1)},
        constrained_layout=True,
    )
    axes[0].barh(y + 0.18, base_values, height=0.34, label="FABDEM brut", color="#7b8fa1")
    axes[0].barh(y - 0.18, hybrid_values, height=0.34, label="MNT hybride", color="#c26a3a")
    axes[0].set_yticks(y, labels)
    axes[0].set_xlabel("Erreur altimétrique (m)")
    axes[0].set_title("Validation croisée spatiale par blocs de 25 km", weight="bold")
    axes[0].grid(axis="x", alpha=0.2)
    axes[0].legend()
    for row, pair in enumerate(zip(base_values, hybrid_values)):
        for offset, value in zip((0.18, -0.18), pair):
            axes[0].text(value + 0.7, row + offset, f"{value:.2f}", va="center", fontsize=9)
    biases = [base["bias_m"], hybrid["bias_m"]]
    colors = ["#7b8fa1", "#c26a3a"]
    axes[1].bar(["FABDEM brut", "MNT hybride"], biases, color=colors)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_ylabel("Biais (m)")
    axes[1].set_title("Biais régional", weight="bold")
    axes[1].grid(axis="y", alpha=0.2)
    for i, value in enumerate(biases):
        axes[1].text(
            i,
            value + (0.09 if value >= 0 else -0.12),
            f"{value:+.2f} m",
            ha="center",
            va="bottom" if value >= 0 else "top",
        )
    fig.suptitle("Apport mesuré du conditionnement par les stations", fontsize=15, weight="bold")
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def method_comparison_figure(final: pd.DataFrame, output: Path) -> None:
    pairs = [
        ("hammer_mgal", "Hammer", "#4c78a8"),
        ("point_mass_mgal", "Masse ponctuelle", "#59a14f"),
        ("vannes_direct_mgal", "Vannes direct (signé)", "#e15759"),
    ]
    reference = final["nagy_prism_mgal"].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.1), constrained_layout=True)
    for ax, (column, label, color) in zip(axes, pairs):
        values = final[column].to_numpy()
        ax.scatter(reference, values, s=8, alpha=0.3, color=color, edgecolors="none")
        lo = min(reference.min(), values.min())
        hi = max(np.percentile(reference, 99.5), np.percentile(values, 99.5))
        ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1, label="égalité")
        corr = np.corrcoef(reference, values)[0, 1]
        mae = np.mean(np.abs(values - reference))
        ax.text(
            0.04,
            0.95,
            f"r = {corr:.3f}\nMAE = {mae:.3f} mGal",
            transform=ax.transAxes,
            va="top",
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.8"},
        )
        ax.set_xlabel("Nagy-prismes (mGal)")
        ax.set_ylabel(f"{label} (mGal)")
        ax.set_title(f"Nagy vs {label}", weight="bold")
        ax.grid(alpha=0.18)
    fig.suptitle(
        "Comparaison station par station des méthodes à 200 km",
        fontsize=15,
        weight="bold",
    )
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def convergence_figure(data: pd.DataFrame, convergence: pd.DataFrame, output: Path) -> None:
    labels = {
        "hammer_mgal": "Hammer",
        "nagy_prism_mgal": "Nagy-prismes",
        "point_mass_mgal": "Masse ponctuelle",
        "vannes_direct_mgal": "Vannes direct",
    }
    colors = ["#4c78a8", "#f28e2b", "#59a14f", "#e15759"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.3), constrained_layout=True)
    for (column, label), color in zip(labels.items(), colors):
        grouped = data.groupby("radius_km")[column]
        radii = np.asarray(sorted(data.radius_km.unique()))
        means = grouped.mean().reindex(radii).to_numpy()
        p95 = grouped.quantile(0.95).reindex(radii).to_numpy()
        axes[0].plot(radii, means, marker="o", color=color, label=f"{label} — moyenne")
        axes[0].plot(
            radii,
            p95,
            marker=".",
            linestyle="--",
            color=color,
            alpha=0.75,
            label=f"{label} — P95",
        )
    axes[0].set_xscale("log")
    axes[0].set_xticks([20, 50, 100, 200], ["20", "50", "100", "200"])
    axes[0].set_xlabel("Rayon externe (km)")
    axes[0].set_ylabel("Correction / effet (mGal)")
    axes[0].set_title("Valeurs cumulées", weight="bold")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8, ncol=2)
    intervals = ["20–50", "50–100", "100–200"]
    for label, color in zip(labels.values(), colors):
        group = convergence[convergence.method.str.startswith(label.split("-")[0].split()[0])]
        if group.empty:
            lookup = {
                "Hammer": "Hammer-style",
                "Nagy-prismes": "Nagy prism",
                "Masse ponctuelle": "Point mass",
                "Vannes direct": "Vannes direct",
            }
            group = convergence[convergence.method == lookup[label]]
        axes[1].plot(
            intervals,
            group.p95_absolute_increment_mgal.to_numpy(),
            marker="o",
            color=color,
            label=label,
        )
    axes[1].axhline(0.01, color="black", linestyle="--", linewidth=1.2, label="seuil 0,01 mGal")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Couronne ajoutée (km)")
    axes[1].set_ylabel("P95 de |l’incrément| (mGal)")
    axes[1].set_title("Test de convergence radiale", weight="bold")
    axes[1].grid(alpha=0.2, which="both")
    axes[1].legend(fontsize=8)
    fig.suptitle("Influence du rayon d’intégration", fontsize=15, weight="bold")
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def qc_map_figure(
    hybrid: np.ndarray,
    final: pd.DataFrame,
    excluded: list[dict],
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 8), constrained_layout=True)
    ax.imshow(relief_rgb(hybrid), extent=EXTENT, origin="upper", interpolation="bilinear")
    ax.scatter(
        final.longitude_deg,
        final.latitude_deg,
        s=5,
        color="white",
        alpha=0.62,
        label="stations utilisées dans les calculs",
    )
    ex_lon = [item["longitude"] for item in excluded]
    ex_lat = [item["latitude"] for item in excluded]
    ax.scatter(
        ex_lon,
        ex_lat,
        s=52,
        marker="x",
        linewidths=1.8,
        color="#d7191c",
        label="24 contrôles exclus du conditionnement",
    )
    ax.set_xlim(8.35, 11.15)
    ax.set_ylim(3.35, 6.65)
    decorate_map(ax, "Contrôle qualité des altitudes de stations")
    ax.legend(loc="lower right", framealpha=0.9)
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hybrid-dir", type=Path, required=True)
    parser.add_argument("--fabdem-dir", type=Path, required=True)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("outputs/reference/terrain_corrections.csv"),
    )
    parser.add_argument("--surfer", type=Path, default=Path("data/reference/out1_surfer.dat"))
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("data/processed/hybrid_dem/hybrid_dem_report.json"),
    )
    parser.add_argument(
        "--excluded",
        type=Path,
        default=Path("data/processed/hybrid_dem/excluded_controls.json"),
    )
    parser.add_argument(
        "--boundary",
        type=Path,
        default=Path("data/reference/cameroon_boundary.geojson"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/reference/figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    hybrid = mosaic(args.hybrid_dir, "*HYBRID*.tif")
    fabdem = mosaic(args.fabdem_dir, "*FABDEM*.tif")
    data = pd.read_csv(args.results)
    final = data[data.radius_km == data.radius_km.max()].copy()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    excluded = json.loads(args.excluded.read_text(encoding="utf-8"))
    boundary = load_boundary(args.boundary)

    study_area_figure(hybrid, final, boundary, args.output_dir / "study_area_hybrid_dem.jpg")
    dem_comparison_figure(
        fabdem,
        hybrid,
        args.surfer,
        final,
        args.output_dir / "dem_spatial_comparison.jpg",
    )
    validation_figure(report, args.output_dir / "hybrid_dem_validation.jpg")
    method_comparison_figure(final, args.output_dir / "method_pairwise_comparison_200km.jpg")
    convergence = pd.read_csv(args.results.parent / "convergence_statistics.csv")
    convergence_figure(data, convergence, args.output_dir / "radial_convergence_detailed.jpg")
    qc_map_figure(hybrid, final, excluded, args.output_dir / "station_altitude_qc_map.jpg")
    print(f"Created six figures in {args.output_dir}")


if __name__ == "__main__":
    main()
