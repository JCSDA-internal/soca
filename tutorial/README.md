# A tutorial for `soca`
After successfully starting an instance of a singularity container and
compiling `soca`, a `tutorial` directory will be installed under
`<path to build>/soca/tutorial`. It contains all the necessary resources
required to run a low resolution MOM6 forecast model and the `soca`
tutorial applications.

The objective of the tutorial is to provide some examples of the basic steps
required to run some of the `soca` data assimilation applications.
The scenario of the tutorial is as follows,
generate synthetic observations from a perturbed MOM6 forecast and subsequently
assimilate the synthetic observations using the unperturbed background as
initial conditions.

## Content of `tutorial`

The content of `tutorial` is as follows:
```console
├── CMakeLists.txt                       # configure the tutorial       
├── config                               # static yaml configurations
│   ├── 3dvar.yaml                       # 3DVAR configuration example
│   ├── gridgen.yaml                     # `soca` grid generation
│   ├── pert_ic.yaml                     # B-matrix randomization
│   ├── staticb.yaml                     # horizontal convolution
│   └── synthetic_obs.yaml               # generate synthetic obs
├── README.md
├── tutorial_3dvarfgat.sh                # 3DVAR FGAT driver for multiple cycles
├── tutorial_3dvar.sh                    # 3DVAR driver, analysis only
├── tutorial_bump_op.sh                  # initialize the horizontal convolution
├── tutorial_gridgen.sh                  # generate the grid and save it to disk
├── tutorial_make_observations.sh        # generate synthetic observations
├── tutorial_perturb_ic.sh               # perturb IC by randomizing the B-matrix
├── tutorial_plot.py                     # observation and state space plotting tools
├── tutorial_synthetic_observations.py   # generate random locations for synthetic observations
└── tutorial_tools.sh                    # generate dynamic yaml configurations
```

Three of the drivers will need to be run first to initialize all the necessary
resources to run the variational data assimilation examples.


## Running the tutorial applications
`cd` into the tutorial directory where `soca` was built,
``` console
cd <path to build>/soca/tutorial
```

#### Create the data assimilation environment
`tutorial_gridgen.sh`, `tutorial_bump_op.sh`, `tutorial_perturb_ic.sh`
and `tutorial_make_observations.sh`
need to be run in sequence as follows:

``` console
./tutorial_gridgen.sh
./tutorial_bump_op.sh
./tutorial_perturb_ic.sh
./tutorial_make_observations.sh
```

`tutorial_gridgen.sh` will create a file `soca_gridspec.nc` under
`./static/` that will be needed by all the other
`soca` applications.

`tutorial_bump_op` initializes a correlation operator using the NICAS method
in the BUMP package of the SABER repository. One file for each number of processors on the machine will be generated at `./static/bump/`

`tutorial_perturb_ic` randomizes the B-matric defined in the `yaml` configuration
files `./config/staticb.yaml` and `pert_ic.yaml`. It generates a single
perturbed state and saves it as a `MOM6` restart file at `./bkg_pert/MOM.res.nc`.  

`tutorial_make_observations.sh` generates synthetic observations by driving the
`MOM6-solo` model using the `soca_hofx.x` executable. `ioda` observation files are
created at the time and locations specified
in `tutorial_synthetic_observations.py`. This application uses generic
observation operators from the UFO repository. A `obs` directory is created
and should contain the following `ioda` observation files:

```console
├── adt.nc4          # absolute dynamic topography
├── insitu.S.nc4     # salinity profiles
├── insitu.T.nc4     # insitu temperature profiles
├── sss.nc4          # sea surface salinity
└── sst.nc4          # sea surface temperature
```

#### 3DVAR example
To run the 3DVAR example, execute the `tutorial_3dvar.sh` script:
```console
./tutorial_3dvar.sh
```
this script perform a 3D variational minimization using the `soca_var.x`
executable for a 24 hour window using the observations generated above.
The executable takes `config/3dvar.yaml` yaml configuration file as an
argument.
A few figures of surface increments are plotted at the end of the script after
the 3DVAR step is done. The figures are `./scratch_3dvar/incr.ssh.png`,
`./scratch_3dvar/incr.sss.png` and `./scratch_3dvar/incr.sst.png` and represent
increments for sea surface height, sea surface salinity and sea surface temperature,
respectively.

#### 3DVAR FGAT with external model
The 3DVARFGAT tutorial is used to show an example of configuration of a
data assimilation experiment cycling through 3 days.
The data assimilation window is 24 hours and the synthetic observations assimilated are
sea surface temperature, sea surface salinity, insitu temperature and salinity and
absolute dynamic topography.

```console
./soca_3dvarfgat.sh
```

Similarly to the `3DVAR` example, figures of surface increments for outer iterations 1 and 2
can be found in `./scratch_3dvarfgat/incr.[1-2].ssh.png`,
`./scratch_3dvarfgat/incr.[1-2].sss.png` and `./scratch_3dvarfgat/incr.[1-2].sst.png`.

Statistics of global mean absolute error of each observation space assimilated can be found
in `./scratch_3dvarfgat/*global_mae.png`.
