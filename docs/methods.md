# Terrain-effect methods

## Common assumptions

All methods use a homogeneous density of 2,670 kg/m3, station elevation from the
original workbook, WGS 84 horizontal coordinates and the same hybrid/FABDEM
terrain cells. Results are in mGal. The default domain reaches 200 km.

## Nagy rectangular prisms

Every DEM cell is a right rectangular vertical prism between the station level
and mean cell elevation. A hill uses positive rock density above the station
plane. A valley is represented as missing rock below that plane, i.e. negative
density contrast. Both geometries yield a positive classical terrain correction.

The exact vertical attraction is evaluated by the eight-corner closed expression
of Nagy (1966). Terms with a finite zero limit are set explicitly to zero. Unit
tests compare a small distant prism against the point-mass asymptote and verify
the sign symmetry above and below the observer.

## Hammer-style compartments

Cells are grouped by logarithmic annulus and azimuth sector. The area-weighted
mean compartment height is inserted into the sector formula:

```text
TC = G rho theta [r2-r1 + sqrt(r1^2+dh^2) - sqrt(r2^2+dh^2)]
```

The historic Hammer chart does not extend to 200 km. Results beyond its original
zones are therefore called *Hammer-style*, with configurable logarithmic rings,
rather than presented as a literal lookup from the 1939 chart.

## Point masses

Each residual prism is collapsed to a mass at its centroid. This approximation
is computationally simple and converges to the prism solution as horizontal
distance becomes large relative to cell dimensions and height. It is retained as
a diagnostic; it is not selected as the reference terrain correction.

## Vannes spherical residual direct effect

`DNEWTC` is transcribed term by term from Appendix A.7.2 of the supplied Vannes
document. At every cell, the kernel for source height is differenced from the
kernel for station height and integrated with the cell solid angle on a sphere of
radius 6,371 km. Negative elevations are set to zero, matching the source code.

This is a signed residual topographic direct effect that also contains the
condensed-mass term in Vannes' formulation. It is not relabelled as the always
positive classical terrain correction. Comparisons must therefore report both
value and definition.

## Multiresolution cells

Aggregation uses arithmetic mean elevation over complete source-pixel blocks.
Since cell horizontal area scales with the number of source pixels, this
preserves the volume represented by each block. Successive radial annuli use 1,
3, 15 and 60 arc-second cells. There is neither overlap nor gap between tiers.

