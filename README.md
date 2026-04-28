JEDI encapsulation of MOM6

(C) Copyright 2017-2024 UCAR.

This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## Building
For more details about JEDI, including installation see: https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/

1. If building the MOM6 SOCA component: unlike other JEDI projects, this project contains its own bundle, it can be built as follows
```
git clone https://github.com/JCSDA/soca.git
mkdir -p soca_build
cd soca_build
ecbuild ../soca/bundle
cd soca
make -j 4
```

To generate doxygen documentation for the Fortran parts of the code, use the `-DENABLE_SOCA_DOC=ON` flag when running `ecbuild`. Documentation will be generated
in the `soca/docs/html` path of the build directory.

See the [JEDI Documentation](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/) for additional details on how to setup, build, and test JEDI projects.

## Documentation

- [docs/postproc_application.md](docs/postproc_application.md): YAML setup for `soca_anpproc.x`
	ensemble post-processing workflow.
- [src/soca/VariableChange/Soca2Cice/README.md](src/soca/VariableChange/Soca2Cice/README.md):
	sea-ice analysis post-processing details for CICE.
- [src/soca/Utils/incrqc/README.md](src/soca/Utils/incrqc/README.md): increment quality-control
	algorithm and configuration details.
