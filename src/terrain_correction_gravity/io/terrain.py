"""Windowed access to 1-degree GeoTIFF tiles with mass-preserving aggregation."""

from __future__ import annotations

import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
from PIL import Image

from ..constants import DEG_TO_RAD, EARTH_RADIUS


@dataclass(frozen=True)
class TerrainCells:
    longitude_deg: np.ndarray
    latitude_deg: np.ndarray
    height_m: np.ndarray
    east_m: np.ndarray
    north_m: np.ndarray
    width_east_m: np.ndarray
    width_north_m: np.ndarray
    area_m2: np.ndarray
    distance_m: np.ndarray
    azimuth_rad: np.ndarray
    solid_angle_sr: np.ndarray

    def __len__(self) -> int:
        return int(self.height_m.size)


def _tile_name(lat0: int, lon0: int, hybrid: bool) -> str:
    if hybrid:
        return f"N{lat0:02d}E{lon0:03d}_HYBRID_FABDEM_STATIONS_30m.tif"
    return f"N{lat0:02d}E{lon0:03d}_FABDEM_V1-2.tif"


class TiledTerrain:
    """Select hybrid tiles in the station domain and FABDEM outside it."""

    def __init__(
        self,
        hybrid_dir: str | Path,
        fabdem_dir: str | Path,
        max_cache_items: int = 8,
    ) -> None:
        self.hybrid_dir = Path(hybrid_dir)
        self.fabdem_dir = Path(fabdem_dir)
        self.max_cache_items = max_cache_items
        self._cache: OrderedDict[tuple[str, int], np.ndarray] = OrderedDict()

    def tile_path(self, lat0: int, lon0: int) -> Path | None:
        hybrid = self.hybrid_dir / _tile_name(lat0, lon0, hybrid=True)
        if hybrid.exists():
            return hybrid
        fabdem = self.fabdem_dir / _tile_name(lat0, lon0, hybrid=False)
        if fabdem.exists():
            return fabdem
        # The acquisition manifest identifies seven absent tiles as wholly
        # oceanic. They are represented at mean sea level. A coverage audit
        # must be run before calculation so a missing land tile is not hidden.
        return None

    def _array(self, path: Path | None, factor: int) -> np.ndarray:
        key = ("<ocean>" if path is None else str(path), factor)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        if path is None:
            source = np.zeros((3600 // factor, 3600 // factor), dtype=np.float32)
            self._cache[key] = source
            return source
        target_size = 3600 // factor
        try:
            with Image.open(path) as image:
                # Cloud-optimized FABDEM files include internal overviews. For
                # aggregated far-field cells, reading a compatible overview is
                # both faster and more robust than decoding the full 30 m frame.
                selected_frame = 0
                selected_size = 3600
                if factor > 1:
                    candidates: list[tuple[int, int]] = []
                    for frame in range(image.n_frames):
                        image.seek(frame)
                        size = image.size[0]
                        if image.size[0] == image.size[1] and size >= target_size and size % target_size == 0:
                            candidates.append((size, frame))
                    if candidates:
                        selected_size, selected_frame = min(candidates)
                image.seek(selected_frame)
                source = np.asarray(image, dtype=np.float32)
        except OSError as error:
            raise OSError(f"Unable to decode terrain tile {path}") from error
        if source.shape != (selected_size, selected_size):
            raise ValueError(f"Unexpected raster shape {source.shape} for {path}")
        source = np.where(np.isfinite(source), source, 0.0)
        aggregate = selected_size // target_size
        if aggregate > 1:
            source = source.reshape(target_size, aggregate, target_size, aggregate).mean(
                axis=(1, 3), dtype=np.float64
            ).astype(np.float32)
        self._cache[key] = source
        # Aggregated overviews are small and reused by nearly every station;
        # retain all of them. Limit only native 30 m frames, which are about
        # 52 MB each.
        native_keys = [cache_key for cache_key in self._cache if cache_key[1] == 1]
        while len(native_keys) > self.max_cache_items:
            oldest_native = native_keys.pop(0)
            del self._cache[oldest_native]
        return source

    def iter_cells(
        self,
        station_lon_deg: float,
        station_lat_deg: float,
        inner_radius_m: float,
        outer_radius_m: float,
        arc_seconds: int,
    ) -> Iterator[TerrainCells]:
        """Yield cells in an annulus; every source pixel contributes through its mean height."""
        if arc_seconds < 1 or 3600 % arc_seconds:
            raise ValueError("arc_seconds must be a positive divisor of 3600")
        if not (0 <= inner_radius_m < outer_radius_m):
            raise ValueError("Require 0 <= inner_radius_m < outer_radius_m")
        factor = arc_seconds
        lat_margin = outer_radius_m / 110_574.0
        lon_margin = outer_radius_m / (
            111_320.0 * max(math.cos(station_lat_deg * DEG_TO_RAD), 0.01)
        )
        lat_min = station_lat_deg - lat_margin
        lat_max = station_lat_deg + lat_margin
        lon_min = station_lon_deg - lon_margin
        lon_max = station_lon_deg + lon_margin

        for lat0 in range(math.floor(lat_min), math.floor(lat_max) + 1):
            for lon0 in range(math.floor(lon_min), math.floor(lon_max) + 1):
                path = self.tile_path(lat0, lon0)
                raster = self._array(path, factor)
                n = raster.shape[0]
                step_deg = 1.0 / n
                col0 = max(0, int(math.floor((lon_min - lon0) / step_deg)))
                col1 = min(n, int(math.ceil((lon_max - lon0) / step_deg)))
                row0 = max(0, int(math.floor((lat0 + 1.0 - lat_max) / step_deg)))
                row1 = min(n, int(math.ceil((lat0 + 1.0 - lat_min) / step_deg)))
                if row0 >= row1 or col0 >= col1:
                    continue

                lons = lon0 + (np.arange(col0, col1) + 0.5) * step_deg
                lats = lat0 + 1.0 - (np.arange(row0, row1) + 0.5) * step_deg
                lon_grid, lat_grid = np.meshgrid(lons, lats)
                height = raster[row0:row1, col0:col1]

                dlat_rad = (lat_grid - station_lat_deg) * DEG_TO_RAD
                dlon_rad = (lon_grid - station_lon_deg) * DEG_TO_RAD
                mean_lat_rad = 0.5 * (lat_grid + station_lat_deg) * DEG_TO_RAD
                north = EARTH_RADIUS * dlat_rad
                east = EARTH_RADIUS * np.cos(mean_lat_rad) * dlon_rad
                distance = np.hypot(east, north)
                mask = (distance >= inner_radius_m) & (distance < outer_radius_m)
                if not np.any(mask):
                    continue

                lat_rad = lat_grid * DEG_TO_RAD
                step_rad = step_deg * DEG_TO_RAD
                width_north = np.full_like(lat_grid, EARTH_RADIUS * step_rad)
                width_east = EARTH_RADIUS * np.cos(lat_rad) * step_rad
                solid_angle = step_rad * (
                    np.sin(lat_rad + 0.5 * step_rad) - np.sin(lat_rad - 0.5 * step_rad)
                )
                area = EARTH_RADIUS**2 * solid_angle
                yield TerrainCells(
                    longitude_deg=lon_grid[mask].astype(np.float64),
                    latitude_deg=lat_grid[mask].astype(np.float64),
                    height_m=height[mask].astype(np.float64),
                    east_m=east[mask].astype(np.float64),
                    north_m=north[mask].astype(np.float64),
                    width_east_m=width_east[mask].astype(np.float64),
                    width_north_m=width_north[mask].astype(np.float64),
                    area_m2=area[mask].astype(np.float64),
                    distance_m=distance[mask].astype(np.float64),
                    azimuth_rad=np.mod(np.arctan2(east[mask], north[mask]), 2.0 * np.pi),
                    solid_angle_sr=solid_angle[mask].astype(np.float64),
                )
