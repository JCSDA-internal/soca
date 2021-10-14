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

# Intial date of the experiment
date_begin="17 April 2018 00:00:00"

for cycle in {1..3}
do
    echo "*****************************************************"
    echo "Starting cycle for $date_begin"
    echo "*****************************************************"

    # Run the forecast model
    # ----------------------
    echo "-- Forecast"
    mkdir -p fcst                                        # where the output of the forecast will go
    ymdh=$(date '+%Y %m %d %H' -d "$date_begin")         # change date in MOM6 namelist
    ../input.nml.sh $ymdh                                #
    ../forecast.yaml.sh "$date_begin"                    # prepare the yaml configuration for a 24 hour forecast
    mpirun -np 2 ../bin/soca_forecast.x ./forecast.yaml >> forecast.log 2>&1 # run a 24 hour forecast

    # 3Dvar FGAT with a pseudo model
    # ------------------------------
    # This implementation of FGAT does not advance the MOM6 model, instead, backgrounds
    # from the forecast step above are read in lieu of time-steping the physics
    echo "-- Start 3DVAR FGAT:"
    mkdir -p obs_out                                 # observation space directory output
    mkdir -p 3dvar_out                               # DA directory output
    ln -sf ../static/bump .                          # link to the NICAS correlation operator

    for outer_iter in {1..2}
    do
        echo "    - Outer iteration $outer_iter"
        echo "                o 3DVAR FGAT with pseudo model"
        mkdir -p ./obs_out/outeriter_$outer_iter
        ../3dvarfgat.yaml.sh "$date_begin" $outer_iter                  # prepare the yaml configuration
        mpirun -np 2 ../bin/soca_var.x ./3dvarfgat.yaml  >> var.log 2>&1 # run 3dvar fgat for a 24 hour cycle

        # Checkpoint analysis (write a MOM6 restart)
        echo "                o Checkpoint analysis"
        ../checkpoint.yaml.sh "$date_begin"                # prepare the yaml configuration
        mpirun -np 2 ../bin/soca_checkpoint_model.x ./checkpoint.yaml  > checkpoint.log 2>&1

        # Re-forecast starting from the analysis as initial conditions
        echo "                o Re-forecast from analysis"
        mv ./RESTART/MOM.res.nc ./INPUT/MOM.res.nc     # swap initial conditions with 3DVAR analysis
        ../forecast.yaml.sh "$date_begin"              # prepare the yaml configuration for a 24 hour forecast
        mpirun -np 2 ../bin/soca_forecast.x ./forecast.yaml >> reforecast.log 2>&1 # run a 24 hour re-forecast
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
echo "Diagnosis of the 3DVAR FGAT"
echo "*****************************************************"
# Mean absolute error statistics of obs-background
# the below script will output png figures inside of the `scratch` directory.
python ../soca_plot.py
