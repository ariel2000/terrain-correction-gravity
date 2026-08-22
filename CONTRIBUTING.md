# Contributing

1. Do not edit `data/raw/original_stations.xlsx`.
2. Record every data exclusion with source row and reason.
3. Keep physical constants and units explicit.
4. Add an analytic, asymptotic or independently computed test for every kernel
   change.
5. Run `pytest` and `ruff check src tests` before proposing a change.
6. Never commit large third-party DEM tiles to ordinary Git history.

