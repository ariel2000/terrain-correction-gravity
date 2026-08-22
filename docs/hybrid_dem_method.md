# Construction and justification of the hybrid DEM

## 1. Rationale

`out1_surfer.dat` is a 598 x 598 regular grid derived by kriging station
elevations. Bilinear re-extraction at the 1,991 stations gives a high correlation
and an RMSE near 10.6 m for positive elevations, but the same observations were
used to build the grid. That diagnostic is therefore circular. With an effective
spacing of approximately 464 x 556 m, the surface also cannot preserve the
30 m-scale morphology needed in the near field of a terrain correction.

FABDEM v1.2 is used as the base because it derives from Copernicus WorldDEM-30
and removes estimated building and forest-height contributions. This matters in
tropical Cameroon, where a digital surface model can follow canopy rather than
bare earth. Copernicus GLO-30 is retained as an independent control, not averaged
into the final product.

## 2. Inputs and coordinate handling

- 1,991 original rows in `original_stations.xlsx`;
- 1,979 unique coordinate pairs;
- FABDEM v1.2, 1 arc-second, WGS 84;
- Copernicus GLO-30, 1 arc-second, WGS 84;
- calculations of station distances use a local metric approximation at the
  mean latitude of the controls.

The vertical datum of the historic station elevations is not documented in the
source workbook. Consequently, the correction is deliberately smooth and
bounded; it is not an exact enforcement of station heights at DEM pixels.

## 3. Control consolidation

At repeated coordinate pairs, the median station elevation is used to define one
conditioning control. The original rows are preserved. A control is excluded if:

1. the station elevation is 0 m while FABDEM exceeds 20 m;
2. `|station - FABDEM| > 150 m`;
3. repeated coordinates have an elevation spread greater than 5 m.

This produces 1,955 accepted and 24 excluded controls. The thresholds are
explicit modelling choices, recorded in the JSON report, rather than hidden data
edits.

## 4. Residual model

For accepted control `i`:

```text
r_i = h_station,i - h_FABDEM(x_i, y_i)
```

Residuals are interpolated by local inverse-distance weighting on the nearest 32
controls:

```text
w_i(q) = 1 / (d_i(q)^1.5 + 5^1.5)
r_hat(q) = sum_i w_i(q) clip(r_i, -40, +40) / sum_i w_i(q)
```

Distances and softening are in kilometres. Individual influence is limited to
plus or minus 40 m so an isolated datum cannot erase satellite morphology. The
correction is estimated on a 15 arc-second support grid and bicubically upsampled
to 1 arc-second. Outside the convex hull of accepted stations its amplitude is
tapered as:

```text
w(q) = exp[-(d_nearest(q) / 20 km)^2]
```

The final model is:

```text
h_hybrid(q) = h_FABDEM(q) + w(q) r_hat(q)
```

FABDEM land pixels are constrained to non-negative final heights. Original
FABDEM water/sea-level pixels are retained.

## 5. Model selection and validation

Random point-wise cross-validation is spatially optimistic for correlated
topography. The selected validation therefore divides controls into 25 km blocks,
assigns 119 spatial blocks to five folds, and withholds entire blocks.

| Metric | FABDEM baseline | Hybrid prediction |
|---|---:|---:|
| Bias (m) | +1.6733 | -0.0520 |
| MAE (m) | 15.6679 | 15.4409 |
| RMSE (m) | 26.0180 | 25.3602 |
| Median error (m) | +0.5172 | -0.2746 |
| P95 absolute error (m) | 58.5083 | 57.2936 |

Training-control diagnostics are reported separately and must not replace the
block-validation metrics.

## 6. Tile continuity

The correction is computed with a common padded border around each 1-degree
tile before bicubic resampling. Across terrestrial seams, the mean correction
difference is 0.018 m and 99% of edge pixels differ by less than 0.104 m.

## 7. Role in the gravity calculation

The hybrid DEM represents surrounding terrain only. The altitude in the Excel
workbook remains the observation height in all kernels. Within the twelve central
tiles the hybrid DEM is used; beyond them, up to 200 km, unmodified FABDEM is
used because there is no station support for extrapolating the residual field.

The complete executable implementation is `scripts/build_hybrid_dem.py` and the
frozen parameter/result record is `data/processed/hybrid_dem/hybrid_dem_report.json`.

