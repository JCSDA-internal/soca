[![travis_develop](https://travis-ci.com/JCSDA/soca.svg?token=Vu1Csdj6JEdxNw6xXKz8&branch=develop)](http://travis-ci.com/JCSDA/soca)
[![Documentation Status](https://readthedocs.com/projects/jointcenterforsatellitedataassimilation-soca/badge/?version=develop)](https://jointcenterforsatellitedataassimilation-soca.readthedocs-hosted.com/en/develop/?badge=develop)
[![codecov](https://codecov.io/gh/JCSDA/soca/branch/develop/graph/badge.svg?token=uFJ62a68D7)](https://codecov.io/gh/JCSDA/soca)

JEDI encapsulation of MOM6  

(C) Copyright 2017-2020 UCAR.

This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

See the [soca ReadTheDocs](https://jointcenterforsatellitedataassimilation-soca.readthedocs-hosted.com/en/latest/?badge=latest) page for more documentation.

## Building

1. If building as part of the larger SOCA coupled system, see the [soca-bundle](https://github.com/JCSDA/soca-bundle)

2. If building the MOM6 SOCA component: unlike other JEDI projects, this project contains its own bundle, it can be built as follows
```
git clone https://github.com/JCSDA/soca.git
mkdir -p soca_build
cd soca_build
ecbuild ../soca/bundle
cd soca
make -j 4
```

3. If building the same way travis-ci builds the MOM6 SOCA bundle:
```
export MAIN_REPO=soca
export LIB_REPOS="fms gsw mom6 crtm fckit atlas oops saber ioda ufo ioda-converters soca-config"
export BUILD_OPT=""
export BUILD_OPT_CRTM="-DBUILD_CRTM=ON"
export BUILD_OPT_OOPS="-DENABLE_OOPS_TOYMODELS=OFF"
export BUILD_OPT_UFO="-DLOCAL_PATH_TESTFILES_IODA=NONE"
export BUILD_OPT_SOCA="-DSOCA_TESTS_FORC_DEFAULT_TOL=ON -DCRTM_FIX_DIR=../../repo.src/crtm/fix"
export MATCH_REPOS="atlas oops saber ioda ioda-converters ufo soca soca-config"
export LFS_REPOS="crtm"
export REPO_CACHE="/path/to/somewhere/repo.cache"
mkdir -p repo.src
cd repo.src
git clone https://github.com/JCSDA/soca.git
cd ..
./repo.src/soca/.github/travisci/prep.sh
```
Note that ccache needs to be installed or loaded. If testing within singularity, swap ccache with ccache-swig in `./github/travisci/build.sh`.
```
 ./repo.src/soca/.github/travisci/build.sh
```

See the [JEDI Documentation](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/) for additional details on how to setup, build, and test JEDI projects.
