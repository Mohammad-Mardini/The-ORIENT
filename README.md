<div align="justify">

The-ORIENT
----------
The ORIENT (cOsmologically deRIved timE-varyiNg Galactic poTential) is a library for integrating trajectories in time-varying gravitational potentials, derived from Milky Way analogues selected from the Illustris-TNG simulation.


:one: Selection criteria and models
-----------------------------------

Milky Way analogues were selected from the TNG100 simulation box based on halo properties at the last snapshot (representing present day). The selection criteria were detailed in Mardini et al. ([2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...88M/abstract)); briefly:

* M<sub>200</sub> between 0.6 and 2.0×10<sup>12</sup> M<sub>⊙</sub>.
* M<sub>stellar</sub> between 4.5 and 8.3×10<sup>10</sup> M<sub>⊙</sub>.
* triaxiality parameter T < 0.35.
* 40% of star particles have circularity parameter ε > 0.7.
* disk-to-total mass ratio > 0.7

123 halos satisfied the criteria, additional quality control narrowed the number down to 54 models.


:two: Properties of the potentials
----------------------------------

Each potential (corresponding to a Galactic analogue) is modelled as a superposition of a Miyamoto–Nagai disk (corresponding to stellar particles in the cosmological simulation) and an NFW sphere (corresponding to gas and dark matter particles), and has seven time-dependent properties in total: the _θ_ and _φ_ directions of the disk’s normal vector, its mass mass _M_<sub>d</sub>, its length scales _a_ and _b_, the central density _ρ_<sub>0</sub> of the NFW component and its length scale _a_. The first column in each `.dat` file in the data folder is the snapshot number from the Illustris simulation, the second is the cosmic age, and the following are the aforementioned model parameters. All quantities are in the {kiloparsec, solar mass, gigayear} unit system.

The fit procedure is described in detail in appendix A of Mardini et al. ([2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...88M/abstract)).


:three: Installation on fresh Mac
---------------------------------
- Get Xcode [https://apps.apple.com/us/app/xcode/id497799835?mt=12](https://apps.apple.com/us/app/xcode/id497799835?mt=12)
- Get Anaconda [https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)
- Use Anaconda to install the following packages:
  - eigen    (you can install it by `conda install -c conda-forge eigen`) 
  - boost    (you can install it by `conda install -c conda-forge boost`)
  - gsl      (you can install it by `conda install -c conda-forge gsl`)
  - pybind11 (you can install it by `conda install -c conda-forge pybind11`)
- Download and install The-ORIENt:

        git clone https://github.com/Mohammad-Mardini/The-ORIENT.git
        cd The-ORIENT/
        export CPATH=$CONDA_PREFIX/include:$CONDA_PREFIX/include/eigen3:$CPATH
        python setup.py install
        
- Check whether you have installed it properly :sparkles: :camel:
  - `cd`
  - `python`
  - `import orient`
  - `dir(orient.Galaxy)`
 
If everything installed properly then you should see the `get_fit_params` and `get_time_limits` attributes :rocket: 


 
:four: Usage
------------

First one should load a potential model. That should be one of the `.dat` files in the data folder (`$ORIENT_INSTALLATION/data`).


    galaxy = orient.Galaxy('subhalo_498768_parameters_processed.dat')


**Note:** currently on macOS the full path of the data file has to be specified, on Linux the relative path is sufficient. That is due to difficulty compiling with the C++ filesystem library on Mojave.

The integration itself can be performed by creating an `Integrate` object like so (see the constructor’s docstring for more details):

    t_min, t_max = galaxy.get_time_limits() # Or specify manually in Gyr in between these bounds.
    stride = 0.15  # Gyr
    ic = [
        120,  0,   0,      # x,  y,  z  in kpc
        0,    40,  0       # vx, vy, vz in km/s
    ]
    result = orient.Integrate(galaxy, ic, t_min, t_max, stride)

The `Integrate` object has a member `t` that is the time in gigayear of the calculated points along the orbit. The integration results can be accessed with the square brackets operator, or with convenience member functions (always in the {kiloparsec, solar mass, gigayear} unit system). For example, to plot obtain the *x* and *y* coordinates as a function of time:

        x, y = result[:,:2].T

**TODO:** 

(1) the convenience member function should accept a time and return an interpolated position or velocity 

(2) if `kms==True` in the constructor the velocity should be returned in km/s.

</div>
