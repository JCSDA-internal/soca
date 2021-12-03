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

ulimit -s unlimited
ulimit -v unlimited

source ./tutorial_tools.sh

# Create a scratch place and cd into it
create_scratch 'scratch_synobs'

# Prepare soca and MOM6 static files
mom6_soca_static $PWD/..

# Create synthetic obs location
mkdir -p obs_scratch
python ../tutorial_synthetic_observations.py

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Replace the un-perturbed background by the perturbed background
cp ../bkg_pert/MOM.res.nc ./INPUT/

# Create synth obs by sampling the model initialized with the
# perturbed background. This application runs MOM6 in-core.
OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_hofx.x ../config/synthetic_obs.yaml

# Concatenate ioda files that contain the synthetic observations
mkdir -p ../obs
ncrcat -O ./obs_scratch/sst.out_*.nc4 ../obs/sst.nc4
ncrcat -O ./obs_scratch/sss.out_*.nc4 ../obs/sss.nc4
ncrcat -O ./obs_scratch/adt.out_*.nc4 ../obs/adt.nc4
ncrcat -O ./obs_scratch/insitu.S.out_*.nc4 ../obs/insitu.S.nc4
ncrcat -O ./obs_scratch/insitu.T.out_*.nc4 ../obs/insitu.T.nc4
