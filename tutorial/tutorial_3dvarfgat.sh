#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Example of a NOT in-core 3DVAR FGAT
#

set -e

ulimit -s unlimited
ulimit -v unlimited

source ./tutorial_tools.sh

# Create a scratch folder and cd into it
create_scratch 'scratch_3dvarfgat'

# Prepare soca and MOM6 static files
mom6_soca_static $PWD/..

# Static grid and correlation operator
ln -sf ../static/soca_gridspec.nc . # link to the grid
ln -sf ../static/bump .             # link to the NICAS correlation operator

# Intial cycling date
date_begin="17 April 2018 00:00:00"

# Create output directories
mkdir -p fcst  # where the output of the forecast will go
mkdir -p obs_out                                 # observation space directory output

# 3 cycles, 24 hour each
for cycle in {1..3}
do
    echo "*****************************************************"
    echo "Cycle initialized $date_begin"
    echo "*****************************************************"

    # Run the forecast model
    # ----------------------
    echo "-- Forecast"
    ymdh=$(date '+%Y %m %d %H' -d "$date_begin")   # change date in MOM6 namelist
    create_input_nml $ymdh                         # create input.nml valid for $date_begin
    forecast_yaml "$date_begin"                    # prepare the yaml configuration for a 24 hour forecast
    OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_forecast.x ./forecast.yaml >> forecast.log 2>&1 # run a 24 hour forecast

    # 3Dvar FGAT with a pseudo model
    # ------------------------------
    # This implementation of FGAT does not advance the MOM6 model, instead, backgrounds
    # from the forecast step above are read in lieu of time-steping the physics
    echo "-- Start 3DVAR FGAT:"
    mkdir -p 3dvar_out  # create DA directory output
    for outer_iter in {1..2}
    do
        echo "    - Outer iteration $outer_iter"
        echo "            o 3DVAR FGAT with pseudo model"
        mkdir -p ./obs_out/outeriter_$outer_iter
        3dvarfgat_yaml "$date_begin" $outer_iter # prepare the yaml configuration
        OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_var.x ./3dvarfgat.yaml  >> var.log 2>&1 # run 3dvar fgat for a 24 hour cycle

        # Checkpoint analysis (write a MOM6 restart containing the analysis state variable)
        echo "            o Checkpoint analysis"
        checkpoint_yaml "$date_begin"                # prepare the yaml configuration
        OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_checkpoint_model.x ./checkpoint.yaml  > checkpoint.log 2>&1

        # Re-forecast starting from the analysis as initial conditions
        echo "            o Re-forecast from analysis"
        mv ./RESTART/MOM.res.nc ./INPUT/MOM.res.nc  # swap initial conditions with 3DVAR analysis
        forecast_yaml "$date_begin"                 # prepare the yaml configuration for a 24 hour forecast
        OMP_NUM_THREADS=1 mpirun -np 2 ../bin/soca_forecast.x ./forecast.yaml >> reforecast.log 2>&1 # run a 24 hour re-forecast
    done # outer_iter

    # Save some of the output
    # -----------------------
    echo "-- Clean and save some of the DA output"
    dirout=$(date '+%Y%m%d%H' -d "$date_begin")
    mkdir -p $dirout
    mv ./obs_out $dirout
    mv ./3dvar_out/ocn.3dvarfgat.iter*.incr.*.nc ./$dirout
    mv ./3dvar_out/ocn.3dvarfgat.an.*00:00:00Z.nc ./$dirout
    rm -rf 3dvar_out

    # Prepare next cycle
    # ------------------
    echo "-- Prepare the next cycle"
    date_begin=$(date '+%d %B %Y %H:%M:%S' -d "$date_begin  $(date +%Z) + 1 day")
    mv ./RESTART/MOM.res.nc ./INPUT/MOM.res.nc    # final conditions becomes
                                                  # initial conditions for the next cyucle

done # cycle

echo "*****************************************************"
echo "Diagnosis for 3DVAR FGAT"
echo "*****************************************************"
# Global Mean absolute error statistics of obs-background
# the below script will output png figures inside of the `scratch` directory.
python -W ignore ../tutorial_plot.py '3dvarfgat'
printf "Figures : \n `ls *.png`\n"
