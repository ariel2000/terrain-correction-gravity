# Limitations and open uncertainties

1. **Vertical datum:** the workbook does not identify whether station heights are
   orthometric, ellipsoidal or tied to another datum. This is the principal
   unresolved input uncertainty.
2. **Station duplicates:** 12 repeated coordinate records exist. They remain
   separate gravity observations; inconsistent heights are excluded only from
   DEM conditioning.
3. **Zero elevations:** 147 entries equal 0 m. Most coincide with a near-sea-level
   Surfer surface, but at least one is locally inconsistent and is flagged.
4. **Density:** 2,670 kg/m3 is spatially uniform. No geological density model is
   available.
5. **Ocean:** bathymetry and seawater density contrast are not modelled. Missing
   wholly oceanic FABDEM tiles are represented at mean sea level.
6. **Earth geometry:** Vannes is spherical. Prism, Hammer and point-mass methods
   use a local tangent-plane geometry within each station calculation.
7. **Hybrid validation:** the RMSE reduction is modest. The hybrid product should
   be described as bias-conditioned FABDEM, not as a field-validated national DEM.
8. **Surfer comparison:** low station residuals for `out1.dat` are not independent
   validation because the same station heights generated the grid.

