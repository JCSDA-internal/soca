# SOCA 1.0.0 RELEASE NOTES

The Joint Effort for Data assimilation Integration (JEDI) is a development project led by the Joint Center for Satellite Data Assimilation (JCSDA). The purpose of JEDI is to provide a comprehensive generic data assimilation system that will accelerate data assimilation developments and the uptake of new observation platforms for the community.

SOCA 1.0.0 is the initial release of the implementation of JEDI for GFDL's Modular Ocean Model version 6 (MOM6).

## SOURCE CODE

All of the SOCA 1.0.0 software dependencies are available from https://github.com/JCSDA

* The generic data assimilation components **OOPS**, **UFO**, **SABER** and **IODA** provide the central data assimilation capabilities. **OOPS** is the heart of JEDI, defining interfaces and data assimilation applications. **UFO** provides generic observation operators. **IODA** provides the in-memory data access and management and IO for observations. **SABER** provides generic interfaces for covariance matrix operations.
* **MOM6** and **FMS** provided by GFDL but containing modifications related to the build system and the exposition of some of the private variables
* The libraries **ECKIT**, **FCKIT** and **ATLAS** provide general utilities used by all other components.

Descriptions of the various components of JEDI are available [here](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/inside/jedi-components/index.html).

## BUILD SYSTEM

Two modes of building and running SOCA are supported with this release:
* Using a [development container](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/learning/tutorials/level2/index.html).
* Using modules maintained on several [HPC platforms](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/using/jedi_environment/modules.html?highlight=modules). These include NOAA's Hera machine, NASA NCCS's Discover machine and Orion at Mississippi State.
* By self installing all of the dependencies using [JEDI-STACK](https://github.com/JCSDA/jedi-stack). This repository includes everything required to build the system, beginning with installation of the source code compilers.

Users can obtain SOCA with:

`git clone -b 1.0.0 https://github.com/jcsda/soca`

For convenience, the SOCA repository includes ecbuild directives that will clone and build all of the SOCA dependencies.


## DOCUMENTATION AND SUPPORT

Documentation for SOCA can be found at [jedi-docs](https://jedi-docs.jcsda.org). Users are encouraged to explore the documentation for detailed descriptions of the source code, development practices and build systems. Note that users are also encouraged to contribute to the documentation, which can be done by submitting to the JEDI-DOCS repository at https://github.com/JCSDA/jedi-docs. JEDI source code includes DOXYGEN comments and the DOXYGEN generated pages can be found under the documentation for [each component](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/inside/jedi-components/index.html). These pages are useful for learning more about the way the software is structured.

## APPLICATIONS

SOCA is provided with input files for performing a number of applications, from observation simulation to variational data assimilation.
SOCA has been implemented at NOAA-EMC and NASA-GMAO with the UFS and GEOS  MOM6 global or regional, at resolutions ranging from 1 to 1/12th degree. However, the only supported resolution for MOM6 in this
release has 25 vertical levels and a tri-polar horizontal grid discretized in 35 meridional and 72 zonal grid points. This low resolution is meant
to be used for unit testing and running tutorials on workstations or laptops with at least 8GB of RAM.

## POST PROCESSING AND VISUALIZATION

### Analysis fields

SOCA outputs FMS restart files in netCDF format that can be visualized easily with for example [Ncview](http://cirrus.ucsd.edu/~pierce/software/ncview/quick_intro.html) (a netCDF visual browser) or [Panoply](https://www.giss.nasa.gov/tools/panoply/).
Examples of Python utility are provided as part of the SOCA tutorial to visualize horizontal fields.

### Observations

Observation output is available in netCDF following formatting structured defined in IODA. Example of simple python utilities
to read the IODA output is also provided in the tutorial.
