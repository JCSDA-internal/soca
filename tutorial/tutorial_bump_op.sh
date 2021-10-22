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
create_scratch 'scratch_bump'

# Prepare soca and MOM6 static files
mom6_soca_static $PWD/..

# Link to previously generated grid
ln -sf ../static/soca_gridspec.nc .

# Create a NICAS horizontal correlation operator
mkdir -p bump
OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_staticbinit.x ../config/staticb.yaml

# Move bump initialization files
cp -r ./bump ../static/
