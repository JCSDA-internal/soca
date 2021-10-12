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

# Create synthetic obs location
mkdir -p obs_scratch
python ../soca_synthetic_observations.py

# Prepare soca and MOM6 static files
../prep.mom6-soca.static.sh $PWD/..

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Run the forecast model
mkdir -p fcst
mpirun -np 2 ../bin/soca_forecast.x ../config/forecast.yaml

# 3Dvar FGAT with a pseudo model
# (reading files from the forecast step above instead of running MOM6 in core)
mkdir -p obs_out
mkdir -p 3dvar_out
ln -s ../static/bump .
mpirun -np 2 ../bin/soca_var.x ../config/3dvarfgat.yaml
