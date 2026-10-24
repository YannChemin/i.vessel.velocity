# i.vessel.velocity

GRASS GIS addon: speed and heading of vessels from a single Sentinel-2
image, from the displacement of the vessel between bands acquired up to
1 s apart (B02, B08, B03, B04). Input: the four bands and a vector map of
vessel points; output: the points with speed, its standard error, heading,
hull-axis check and wake check.

See `i.vessel.velocity.md` for the method, its precision and its limits.

## Install

Symlinked into `grass-addons/src/imagery`, then:

```sh
make -C ~/dev/grass-addons/src/imagery/i.vessel.velocity MODULE_TOPDIR=$HOME/dev/grass
```

## Tests

```sh
PYTHONPATH=$(grass --config python_path) python3 -m pytest tests
```

The tests build synthetic scenes (moving hulls with wakes) in a temporary
project; no data download.

## Licence

Unlicense (public domain), see `LICENSE`.
