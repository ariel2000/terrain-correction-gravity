# Data dictionary

## Original stations

`data/raw/original_stations.xlsx`, worksheet `Feuil1`:

| Column | Unit | Role |
|---|---|---|
| Longitude | decimal degree east, WGS 84 assumed | observation coordinate |
| Latitude | decimal degree north, WGS 84 assumed | observation coordinate |
| Altitude | m | observation height; never replaced by DEM height |
| A(air libre) | mGal | input free-air anomaly |

Excel row numbers are retained in calculated outputs as `excel_row`.

## Terrain-correction output

| Field | Unit | Meaning |
|---|---|---|
| radius_km | km | external integration radius |
| cell_count | count | cells used after multiresolution selection |
| hammer_mgal | mGal | positive Hammer-style correction |
| nagy_prism_mgal | mGal | positive exact-prism correction |
| point_mass_mgal | mGal | positive centroid approximation |
| vannes_direct_mgal | mGal | signed Vannes spherical residual direct effect |

## Integrity hashes

| File | SHA-256 |
|---|---|
| original_stations.xlsx | `167c8d10864b3099adc95fb340f1322fa164034c4b308c9d0318be28eb6a19e6` |
| out1_surfer.dat | `273622677bc05343f87e38c0141bded3e36c9b87c2b773aec9af63b3ed4e2e2d` |
| Vannes appendix PDF | `a45e1f4cb3e405d01128d3e2e3d45e3a927161db0916356459b7a7b0ffb3ef4c` |

