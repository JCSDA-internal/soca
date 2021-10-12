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

# Create a scratch place
[ -d scratch ] && rm -rf scratch
mkdir scratch
cd scratch

# Prepare soca and MOM6 static files
../prep.mom6-soca.static.sh $PWD/..

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Link to previously generated NICAS horizontal correlation operator
ln -sf ../static/bump .

# Generate a perturbation by randomizing a staic B-matrix
mkdir -p out
mpirun -np 2 ../bin/soca_enspert.x ../config/ensemble_perturbations.yaml

# Move the perturbed restarts
mkdir -p ../bkg_pert
mv ./RESTART/MOM.res.nc ../bkg_pert
