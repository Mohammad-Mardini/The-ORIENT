# The ORIENT

**The ORIENT** (c**O**smologically de**RI**ved tim**E**-varyi**N**g Galactic po**T**ential) is a library for integrating orbits in time-varying gravitational potentials. The potentials are derived from Milky Way analogues selected from the IllustrisTNG cosmological simulation (TNG100).

## Contents

1. [Milky Way analogues](#1-milky-way-analogues)
2. [Potential models](#2-potential-models)
3. [Installation](#3-installation)
4. [Usage](#4-usage)
5. [Controlling the time dependence](#5-controlling-the-time-dependence)
6. [Citation](#6-citation)
7. [Roadmap](#7-roadmap)

---

## 1. Milky Way analogues

Milky Way analogues were selected from the TNG100 box using halo properties at the final snapshot (present day). The full selection is described in [Mardini et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...88M/abstract). In brief:

| Property | Criterion |
|---|---|
| Halo mass, M<sub>200</sub> | 0.6–2.0 × 10<sup>12</sup> M<sub>⊙</sub> |
| Stellar mass, M<sub>★</sub> | 4.5–8.3 × 10<sup>10</sup> M<sub>⊙</sub> |
| Triaxiality parameter, T | < 0.35 |
| Circularity | ≥ 40% of star particles with ε > 0.7 |
| Disk-to-total mass ratio | > 0.7 |

Of the 123 halos that satisfy these criteria, 54 passed additional quality control and are provided as potential models.

## 2. Potential models

Each potential is the superposition of two components:

- a **Miyamoto–Nagai disk**, fit to the stellar particles;
- an **NFW sphere**, fit to the gas and dark matter particles.

Each model has seven time-dependent parameters. Every `.dat` file in the `data/` folder has one row per simulation snapshot, with the following columns:

| Column | Quantity |
|---|---|
| 1 | Illustris snapshot number |
| 2 | Cosmic age |
| 3 | θ, polar angle of the disk normal vector |
| 4 | φ, azimuthal angle of the disk normal vector |
| 5 | M<sub>d</sub>, disk mass |
| 6 | a, disk radial scale length |
| 7 | b, disk vertical scale length |
| 8 | ρ<sub>0</sub>, NFW central density |
| 9 | NFW scale length |

All quantities use the {kpc, M<sub>⊙</sub>, Gyr} unit system. The fitting procedure is described in Appendix A of [Mardini et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...88M/abstract).

## 3. Installation

### Requirements

- Python 3 with NumPy
- A C++ compiler (on macOS, install [Xcode](https://apps.apple.com/us/app/xcode/id497799835) or the Xcode Command Line Tools)
- [Eigen](https://eigen.tuxfamily.org), [Boost](https://www.boost.org), [GSL](https://www.gnu.org/software/gsl/), and [pybind11](https://github.com/pybind/pybind11)

The simplest route is a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment:

```bash
conda install -c conda-forge eigen boost gsl pybind11 numpy
```

### Build and install

```bash
git clone https://github.com/Mohammad-Mardini/The-ORIENT.git
cd The-ORIENT
export CPATH=$CONDA_PREFIX/include:$CONDA_PREFIX/include/eigen3:$CPATH
pip install . --no-build-isolation
```

### Check the installation

From any directory other than the repository:

```bash
python -c "import orient"
```

No output means the installation succeeded.

## 4. Usage

### Load a potential

Choose one of the `.dat` files in `data/`:

```python
import numpy as np
import orient

galaxy = orient.Galaxy('data/subhalo_498768_parameters_processed.dat')
```

> **Note (macOS):** on some macOS versions the full path to the data file is required, because the C++ `filesystem` library could not be used when compiling on Mojave. On Linux a relative path works.

### Integrate an orbit

Create an `Integrate` object (see its docstring for all options):

```python
t_min, t_max = galaxy.get_time_limits()   # or any range within these bounds (Gyr)
stride = -8.0                              # Gyr

ic = [120.0,  0.0, 0.0,    # x,  y,  z  [kpc]
        0.0, 40.0, 0.0]    # vx, vy, vz [km/s]

result = orient.Integrate(galaxy, ic, t_min, t_max, stride)
```

### Access the results

`result.t` holds the times (Gyr) of the computed orbit points. Positions are available as convenience members, and all phase-space coordinates through the `[]` operator. Outputs are in {kpc, M<sub>⊙</sub>, Gyr} units.

```python
t = result.t
x, y, z = result.x, result.y, result.z
```

### Orbital parameters

The library also computes orbital parameters such as eccentricity and maximum height above the disk plane, evaluated at chosen look-back times:

```python
lookback = np.array([-2, -4, -6, -8])     # Gyr

eccentricity = result.e(t_max + lookback)
z_max        = result.zmax(t_max + lookback)
```

### Example

The animation below shows the integrated orbit of the primordial thin-disk star 2MASS J1808−5104 over the past 8 Gyr. See [Mardini et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmp.2594M/abstract) for details. A tutorial combining The ORIENT with Gaia astrometry is included in the repository.

https://user-images.githubusercontent.com/35367221/195906213-2de68c3b-3d70-431f-a21f-8d0fbbf06622.mp4

## 5. Controlling the time dependence

Any parameter can be held fixed at its value at a given time. To turn The ORIENT into a static potential at the present day:

```python
galaxy.freeze_all_parameters(t_max)
```

To freeze only the orientation of the disk:

```python
galaxy.freeze_parameter('theta', t_max)
galaxy.freeze_parameter('phi', t_max)
```

## 6. Citation

If you use The ORIENT in your work, please cite:

- Mardini, M. K., et al. 2020, ApJ, 903, 88 · [ADS](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...88M/abstract)
- Mardini, M. K., et al. 2022, MNRAS · [ADS](https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmp.2594M/abstract)

## 7. Roadmap

- [ ] The convenience member function should accept a time and return an interpolated position or velocity 
- [ ] Return velocities in km/s when `kms=True` is passed to the `Integrate` constructor.

  
