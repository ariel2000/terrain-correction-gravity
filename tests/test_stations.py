from pathlib import Path

from terrain_correction_gravity.io.stations import audit_stations, load_stations

ROOT = Path(__file__).resolve().parents[1]


def test_original_station_inventory():
    stations = load_stations(ROOT / "data/raw/original_stations.xlsx")
    audit = audit_stations(stations)
    assert audit["count"] == 1991
    assert audit["unique_coordinate_count"] == 1979
    assert audit["duplicate_coordinate_count"] == 12
    assert audit["zero_altitude_count"] == 147
    assert audit["extent_wgs84"] == {
        "west": 8.502,
        "east": 10.9933,
        "south": 3.5006,
        "north": 6.5,
    }
