# Results

The complete run contains 7,964 station-radius rows (1,991
stations and 4 radii). All reported method values are finite.

## Distribution at 200 km

| Method | Mean | Median | P95 | Maximum | Negative stations |
|---|---:|---:|---:|---:|---:|
| Hammer-style | 2.006 | 0.864 | 7.463 | 61.328 | 0 |
| Nagy prism | 3.217 | 1.906 | 10.478 | 68.556 | 0 |
| Point mass | 3.045 | 1.532 | 10.738 | 64.931 | 0 |
| Vannes direct | 2.880 | 1.306 | 12.235 | 158.643 | 400 |

## Convergence from 100 to 200 km

| Method | Mean increment | P95 absolute increment | Fraction below 0.01 mGal |
|---|---:|---:|---:|
| Hammer-style | 0.089 | 0.287 | 3.2% |
| Nagy prism | 0.102 | 0.302 | 1.0% |
| Point mass | 0.102 | 0.302 | 1.0% |
| Vannes direct | -0.027 | 0.481 | 2.9% |


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
