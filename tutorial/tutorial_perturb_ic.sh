#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Create a perturbed state by randomizing a static B-matrix
#

set -e

ulimit -s unlimited
ulimit -v unlimited

source ./tutorial_tools.sh

# Create a scratch place and cd into it
create_scratch 'scratch_pertic'

# Prepare soca and MOM6 static files
mom6_soca_static $PWD/..

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Link to previously generated NICAS horizontal correlation operator
ln -sf ../static/bump .

# Generate a perturbation by randomizing a staic B-matrix
mkdir -p out
OMP_NUM_THREADS=1 mpirun ../bin/soca_enspert.x ../config/pert_ic.yaml

# Move the perturbed restarts
mkdir -p ../bkg_pert
mv ./RESTART/MOM.res.nc ../bkg_pert
