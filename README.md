# Minor Planet Painter 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
A Python library for plotting and visualizing minor planets


## Structure
```
./
  data/
  minor_planet_painter/
    common.py
    ...
  scripts/
    plot_sssb_xy.py
    ...
  .gitignored
  README.md
```

0. Obtain datasets from MPC
``` 
# Orbital elements of all minorplanets (MPCORB.DAT) and NEAs (NEAm00.txt) are saved in ./data
wget_MPCORB_NEA.sh
```

1. Spatial distribution of minor bodies
```
# Plot all minor planets
plot_sssbs_xy.py

# Plot all minor planets specifying the input file
plot_asteroid_xy.py --

```

## Installing
```
git clone git@github.com:jinbeniyama/minor-planet-painter.git
```
