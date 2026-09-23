# Minor Planet Painter 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
A Python library for plotting and visualizing minor planets


## Structure
```
./
  data/
  fig/
  minor_planet_painter/
    common.py
    ...
  scripts/
    plot_sssb_xy.py
    ...
  .gitignored
  README.md
```

## Set up
```
pipx install --editable . --force
```

## Make figures
0. Preparation

- MPC
Obtain datasef from MPC. Orbital elements of all minorplanets (`MPCORB.DAT`) and NEAs (`NEAm00.txt`) are saved in `./data`.
``` 
wget_MPCORB_NEA.sh
```

- LCDB
You can download from the [LCDB website](https://minplanobs.org/MPInfo/php/lcdb.php).
Download The latest data release, 2023 Oct, from [here](https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip).
The input file used in this repository is `lc_summary_pub.txt`.


1. Spatial distribution of minor bodies
```
# Plot all minor planets (output figure is shown below)
plot_sssb_xy.py 2025-08-25 --range 60
# Plot all minor planets, only inner system
plot_sssb_xy.py 2025-08-25 --range 6

# Plot all minor planets specifying the input file
plot_sssb_xy.py --MPCORB MPCORB_original.DAT
```

![Spatial distribution of minor bodies](fig/MPCORB_20250825.jpg)


2. Orbital elements of minor bodies
```
# Plot only NEAs (output figure is shown below)
plot_sssb_orbelem.py --onlyNEA

# NEA pairs (in prep)
```

![Orbital elements of minor bodies](fig/e_NEA.jpg)


3. Angular distance of minor bodies
```
# Plot all minor planets (output figure is shown below)
# Not that this is from the Sun, not the Earth.
plot_sssb_angsize.py 2025-08-25 --out angsize_20250825.jpg
```

![Angular size of minor bodies](fig/angsize_20250825.jpg)

4. Sky motion of minor bodies
```
# Plot all minor planets (output figure is shown below)
# It takes a few minutes for only 500 bodies
plot_sssb_skymotion.py --out skymotion_20250825.jpg --Nobj 500
```
<p align="center">
  <img src="/fig/skymotion_20250825.jpg" width="600"/><br>
  <em>Sky motion of minor bodies. The NEA with a semimajor axis of about 1 is (99942) Apophis.</em>
</p>


5. Psid vs. Psyn
```
plot_sssb_Psid_Psyn.py
```

6. Lightcurve database
```
plot_lcdb_animation.py data/LCLIST_PUB_CURRENT_2023Oct/lc_summary_pub.txt
```

## Installing
```
git clone git@github.com:jinbeniyama/minor-planet-painter.git
```
