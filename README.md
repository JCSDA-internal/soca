[![travis_develop](https://travis-ci.com/JCSDA/soca.svg?token=Vu1Csdj6JEdxNw6xXKz8&branch=develop)](http://travis-ci.com/JCSDA/soca)
[![Documentation Status](https://readthedocs.com/projects/jointcenterforsatellitedataassimilation-soca/badge/?version=develop)](https://jointcenterforsatellitedataassimilation-soca.readthedocs-hosted.com/en/develop/?badge=develop)
[![codecov](https://codecov.io/gh/JCSDA/soca/branch/develop/graph/badge.svg?token=uFJ62a68D7)](https://codecov.io/gh/JCSDA/soca)

JEDI encapsulation of MOM6  

(C) Copyright 2017-2020 UCAR.

This software is licensed under the terms of the Apache Licence Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

See the [soca ReadTheDocs](https://jointcenterforsatellitedataassimilation-soca.readthedocs-hosted.com/en/latest/?badge=latest) page for more documentation.

## Building

1. If building as part of the larger SOCA coupled system, see the [soca-bundle](https://github.com/JCSDA/soca-bundle)

2. Otherwise, for just the MOM6 SOCA component: unlike other JEDI projects, this project contains its own bundle, it can be built as follows
```
git clone https://github.com/JCSDA/soca.git
mkdir -p soca_build
cd soca_build
ecbuild ../soca/bundle
cd soca
make -j 4
```

See the [JEDI Documentation](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/) for additional details on how to setup, build, and test JEDI projects.