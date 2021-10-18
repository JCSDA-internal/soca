#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# 3DVAR example
#

set -e

ulimit -s unlimited
ulimit -v unlimited

source ./tutorial_tools.sh

# Create a scratch folder and cd into it
create_scratch 'scratch_3dvar'

# Prepare soca and MOM6 static files
mom6_soca_static $PWD/..

# Static grid and correlation operator
ln -sf ../static/soca_gridspec.nc . # link to the grid
ln -sf ../static/bump .             # link to the NICAS correlation operator

# Create output directories
mkdir -p 3dvar_out
mkdir -p obs_out

# 3DVAR
OMP_NUM_THREADS=1 mpirun ../bin/soca_var.x ../config/3dvar.yaml

# Plot some diagnostics
python -W ignore ../tutorial_plot.py '3dvar'
printf "Figures : \n `ls *.png`\n"
