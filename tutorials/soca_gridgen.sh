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

# Create a scratch place
[ -d scratch ] && rm -rf scratch
mkdir scratch
cd scratch

# Prepare soca and MOM6 static files
../prep.mom6-soca.static.sh $PWD/..

# Generate grid
mpirun -np 2 ../bin/soca_gridgen.x ../config/gridgen.yaml

# Save grid for later use
mkdir -p ../static
mv ./soca_gridspec.nc ../static/
