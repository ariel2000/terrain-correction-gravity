"""Shared-domain calculation pipeline for all terrain-effect methods."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io.terrain import TerrainCells, TiledTerrain
from .methods.hammer import hammer_sector_correction, logarithmic_edges
from .methods.point_mass import point_mass_correction
from .methods.prism import nagy_prism_correction
from .methods.vannes import vannes_adapted_correction, vannes_direct_effect


@dataclass(frozen=True)
class ResolutionTier:
    inner_km: float
    outer_km: float
    arc_seconds: int


def concatenate_cells(parts: list[TerrainCells]) -> TerrainCells:
    if not parts:
        empty = np.array([], dtype=np.float64)
        return TerrainCells(*(empty.copy() for _ in range(11)))
    fields = TerrainCells.__dataclass_fields__
    return TerrainCells(
        **{name: np.concatenate([getattr(part, name) for part in parts]) for name in fields}
    )


def subset_cells(cells: TerrainCells, mask: np.ndarray) -> TerrainCells:
    return TerrainCells(
        **{name: getattr(cells, name)[mask] for name in TerrainCells.__dataclass_fields__}
    )


def collect_cells(
    terrain: TiledTerrain,
    longitude_deg: float,
    latitude_deg: float,
    radius_km: float,
    tiers: list[ResolutionTier],
) -> TerrainCells:
    parts: list[TerrainCells] = []
    maximum_m = radius_km * 1000.0
    for tier in tiers:
        inner_m = tier.inner_km * 1000.0
        outer_m = min(tier.outer_km * 1000.0, maximum_m)
        if outer_m <= inner_m:
            continue
        parts.extend(
            terrain.iter_cells(
                longitude_deg,
                latitude_deg,
                inner_m,
                outer_m,
                tier.arc_seconds,
            )
        )
    return concatenate_cells(parts)


def calculate_station(
    terrain: TiledTerrain,
    longitude_deg: float,
    latitude_deg: float,
    station_height_m: float,
    radius_km: float,
    tiers: list[ResolutionTier],
    density_kg_m3: float = 2670.0,
    hammer_sectors: int = 12,
    hammer_rings_per_decade: int = 8,
) -> dict[str, float | int]:
    cells = collect_cells(terrain, longitude_deg, latitude_deg, radius_km, tiers)
    edges = logarithmic_edges(radius_km * 1000.0, rings_per_decade=hammer_rings_per_decade)
    return {
        "cell_count": len(cells),
        "hammer_mgal": hammer_sector_correction(
            cells, station_height_m, edges, hammer_sectors, density_kg_m3
        ),
        "nagy_prism_mgal": nagy_prism_correction(cells, station_height_m, density_kg_m3),
        "point_mass_mgal": point_mass_correction(cells, station_height_m, density_kg_m3),
        "vannes_direct_mgal": vannes_direct_effect(
            cells,
            longitude_deg,
            latitude_deg,
            station_height_m,
            density_kg_m3,
        ),
        "vannes_adapted_mgal": vannes_adapted_correction(
            cells,
            longitude_deg,
            latitude_deg,
            station_height_m,
            density_kg_m3,
        ),
    }


def calculate_station_radii(
    terrain: TiledTerrain,
    longitude_deg: float,
    latitude_deg: float,
    station_height_m: float,
    radii_km: list[float],
    tiers: list[ResolutionTier],
    density_kg_m3: float = 2670.0,
    hammer_sectors: int = 12,
    hammer_rings_per_decade: int = 8,
) -> list[dict[str, float | int]]:
    """Collect the largest domain once, then evaluate all convergence radii."""
    if not radii_km:
        return []
    maximum = max(radii_km)
    all_cells = collect_cells(terrain, longitude_deg, latitude_deg, maximum, tiers)
    outputs: list[dict[str, float | int]] = []
    for radius in radii_km:
        cells = subset_cells(all_cells, all_cells.distance_m < radius * 1000.0)
        edges = logarithmic_edges(
            radius * 1000.0, rings_per_decade=hammer_rings_per_decade
        )
        outputs.append(
            {
                "cell_count": len(cells),
                "hammer_mgal": hammer_sector_correction(
                    cells, station_height_m, edges, hammer_sectors, density_kg_m3
                ),
                "nagy_prism_mgal": nagy_prism_correction(
                    cells, station_height_m, density_kg_m3
                ),
                "point_mass_mgal": point_mass_correction(
                    cells, station_height_m, density_kg_m3
                ),
                "vannes_direct_mgal": vannes_direct_effect(
                    cells,
                    longitude_deg,
                    latitude_deg,
                    station_height_m,
                    density_kg_m3,
                ),
                "vannes_adapted_mgal": vannes_adapted_correction(
                    cells,
                    longitude_deg,
                    latitude_deg,
                    station_height_m,
                    density_kg_m3,
                ),
            }
        )
    return outputs
