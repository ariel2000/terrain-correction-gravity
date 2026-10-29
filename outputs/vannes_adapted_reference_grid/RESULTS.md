# Controlled SCCV benchmark

The benchmark contains 1,991 gravity stations evaluated on the common
information-equivalent 200 x 200 historical terrain grid.  Every station uses
40,000 terrain cells, a density of 2,670 kg/m3 and the complete rectangular
grid domain.  SCCV uses exact Nagy prisms within 5 km and an order-eight
Gauss-Legendre spherical integration outside that radius.

## Comparison with Nagy prisms

| Indicator | Original Vannes | SCCV |
|---|---:|---:|
| Mean bias (mGal) | -1.866119 | -0.002618 |
| MAE (mGal) | 2.323545 | 0.002619 |
| RMSE (mGal) | 5.098439 | 0.006768 |
| P95 absolute difference (mGal) | 9.448417 | 0.011113 |
| Maximum absolute difference (mGal) | 73.950194 | 0.127825 |
| Pearson correlation | -0.274040 | 0.999999 |
| Negative corrections | 604 | 0 |

For SCCV, 93.97% of stations are within 0.01 mGal of Nagy, 99.75% are
within 0.05 mGal and 99.95% are within 0.10 mGal.

## Interpretation

The result demonstrates numerical and physical coherence under controlled
common conditions.  It is not a fully independent validation because SCCV and
the Nagy reference share the exact prism calculation within 5 km.  Independent
tesseroid or three-dimensional spherical volume integration, sensitivity to
the near-field radius, vertical-datum verification and regional radial
convergence remain required before claiming definitive validation.

The script generates `station_comparison.csv` locally with pseudonymised
station-level calculations.  That file is excluded from the public repository
pending signed data-publication authorisation. `pairwise_summary.csv` contains
the non-identifying reproducible comparison statistics.
