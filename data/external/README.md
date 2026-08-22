# External DEM files

Large third-party GeoTIFF files are not committed to Git.

Expected directories:

```text
data/external/fabdem_v1_2/N03E008_FABDEM_V1-2.tif
data/external/copernicus_glo30/Copernicus_DSM_COG_10_N03_00_E008_00_DEM.tif
data/processed/hybrid_dem/tiles/N03E008_HYBRID_FABDEM_STATIONS_30m.tif
```

Coverage acquired for this study is approximately 6.70-12.80 deg E and
1.69-8.31 deg N, sufficient for a 200 km calculation around the stations:

- 49 land-containing FABDEM v1.2 tiles;
- 49 corresponding Copernicus GLO-30 tiles;
- seven absent tiles in both products, identified as wholly oceanic;
- 12 hybrid tiles N03E008 through N06E010 in the central station domain.

`tile_manifest.csv` records the filename, byte size and SHA-256 digest of all
98 acquired source GeoTIFFs, plus the 14 expected-but-absent ocean entries.

Sources:

- FABDEM v1.2: https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn
- Copernicus DEM GLO-30: https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM

FABDEM attribution and non-commercial conditions must be retained. Copernicus
terms also apply independently of this repository's MIT software license.
