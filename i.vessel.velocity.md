## DESCRIPTION

*i.vessel.velocity* estimates the speed and heading of vessels from a
single multispectral image of a push-broom sensor, by default the
Sentinel-2 MSI. The MSI records each band at a slightly different time:
the green band (B03) 0.527 s, the red band (B04) 1.005 s and the
near-infrared band (B08) 0.264 s after the blue band (B02). The bands of
a Level-2A product are orthorectified for the ground, so a static target
sits at the same place in every band, while a moving vessel shifts by its
velocity times the time lag (Heiselberg 2019). A vessel at 10 knots moves
about 5 m, half a 10 m pixel, between B02 and B04.

For each point of the **input** vector map, usually the barycentres of
detected vessels, the module:

1. reads a small window of the four bands around the point;
2. subtracts the local sea background (median of a ring around the hull)
   and computes the intensity-weighted barycentre of the bright blob in
   each band;
3. fits position against acquisition time by least squares, which gives
   the velocity components, the standard error of the speed from the
   scatter of the four positions, and the heading as a bearing from grid
   north;
4. compares the heading with the hull axis, if **orient_column** is given:
   a vessel under way moves along its long axis;
5. for a moving vessel, compares the mean near-infrared brightness on a
   strip one to three hull lengths behind it with the same strip ahead:
   the wake trails behind, so a larger value behind supports the direction
   of the heading.

The **output** map is a copy of the input with these columns added:

| Column | Content |
|---|---|
| `n_bands` | bands in which the vessel was found (3 or 4 to measure) |
| `snr_min` | lowest peak signal-to-noise ratio over the bands |
| `vx_ms`, `vy_ms` | velocity, east and north, m/s |
| `speed_kn`, `speed_se_kn` | speed and its standard error, knots |
| `heading_deg` | heading, degrees clockwise from grid north |
| `axis_diff_deg` | angle between heading and hull axis, 0-90 degrees |
| `moving` | 1 when the speed exceeds **sigma** standard errors and **min_speed** |
| `wake_behind`, `wake_ahead` | NIR brightness behind and ahead (moving vessels) |
| `kin_status` | `ok`, `few_bands` (found in fewer than 3 bands) or `edge` |

The measurement window holds the hull plus three hull lengths on each
side; points whose window leaves the computational region are flagged
`edge`. The region must match the band grid.

## NOTES

The precision is set by the barycentre of a hull a few pixels long: on
Sentinel-2 vessels of 25-40 m, the standard error of the speed is
typically 1-3 knots, so only vessels above about 5-8 knots are reliably
flagged as moving. Larger hulls are measured better.

The time lags of the MSI vary by up to 2 % across and along the swath
(Binet et al. 2022), which is small next to this precision. The
registration between the bands is not: the Sentinel-2 requirement is
0.3 pixel at the 3-sigma level and ESA measures 0.07-0.22 pixel for the
band pairs (Sentinel-2 Data Quality Report), and 0.2 pixel, 2 m, read as
a displacement over 1 s, is about 4 knots. A systematic apparent velocity found on stationary vessels of
the same orbit and satellite (for example ships at anchor) can be removed
with **bias**. Registration measured on land is not a substitute: the
contrast of vegetation and terrain shadows differs between the visible
and near-infrared bands and shifts their apparent position.

The heading is relative to grid north; the difference from true north
(the meridian convergence of the projection) is a fraction of a degree
near the central meridian of a UTM zone.

Other sensors with a known inter-band delay can be used by giving their
bands and **time_lags**, in the order blue, green, red, near-infrared.

## EXAMPLES

Speed and heading of the vessels detected in one Sentinel-2 scene,
imported with *r.in.sentinel*:

```sh
r.in.sentinel -o collection=sentinel-2-l2a bands=B02,B03,B04,B08 \
  start=2024-10-20 end=2024-10-21 resolution=10 output=s2
g.region raster=s2_20241020_B08

i.vessel.velocity input=vessels blue=s2_20241020_B02 green=s2_20241020_B03 \
  red=s2_20241020_B04 nir=s2_20241020_B08 output=vessels_velocity \
  length_column=length_m orient_column=orient_deg

v.db.select vessels_velocity columns=cat,speed_kn,speed_se_kn,heading_deg,moving \
  where="moving = 1"
```

## REFERENCES

- Heiselberg, H. (2019). Aircraft and ship velocity determination in
  Sentinel-2 multispectral images. *Sensors*, 19(13), 2873.
  [doi:10.3390/s19132873](https://doi.org/10.3390/s19132873)
- Binet, R., Bergsma, E., Poulain, V. (2022). Accurate Sentinel-2
  inter-band time delays. *ISPRS Annals of the Photogrammetry, Remote
  Sensing and Spatial Information Sciences*, V-1-2022, 57-66.
  [doi:10.5194/isprs-annals-V-1-2022-57-2022](https://doi.org/10.5194/isprs-annals-V-1-2022-57-2022)
- ESA, Sentinel-2 MSI technical guide: MSI instrument (inter-band
  temporal offsets).

## SEE ALSO

*[i.landsat.toar](https://grass.osgeo.org/grass-stable/manuals/i.landsat.toar.html),
[r.in.sentinel](https://github.com/YannChemin/r.in.sentinel),
[v.db.select](https://grass.osgeo.org/grass-stable/manuals/v.db.select.html),
[v.to.db](https://grass.osgeo.org/grass-stable/manuals/v.to.db.html)*

## AUTHORS

Yann Chemin
