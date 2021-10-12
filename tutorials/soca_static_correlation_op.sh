#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Make synthetic observations for the soca tutorials.
#

set -e

soca_src=/home/gvernier/sandboxes/SOCA-1.0.0/soca-release
soca_build=/home/gvernier/sandboxes/SOCA-1.0.0/build.soca/

# Create a scratch place
[ -d scratch ] && rm -rf scratch
mkdir scratch
cd scratch

# Prepare soca and MOM6 static files
../prep.mom6-soca.static.sh $PWD/..

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Create a NICAS horizontal correlation operator
mkdir -p bump
mpirun -np 2 ../bin/soca_staticbinit.x ../config/staticb.yaml

# Move bump initialization files
mv ./bump ../static/
