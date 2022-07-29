#!/bin/bash

mkdir -p ./INPUT_small/
cp ./INPUT/hycom1_25.nc ./INPUT_small/
cp ./INPUT/layer_coord25.nc ./INPUT_small/
ncks -O -d nxp,0,144,2 -d nx,0,143,2 -d nyp,0,69,2 -d ny,0,67,2 ./INPUT/ocean_hgrid.nc ./INPUT_small/ocean_hgrid_small.nc
ncks -O -d nx,0,71,2 -d ny,0,33,2 ./INPUT/ocean_topog.nc ./INPUT_small/ocean_topog_small.nc 
ncks -O -d lonq,0,71,2 -d lonh,0,71,2 -d lath,0,33,2 -d latq,0,33,2 ./INPUT/MOM.res.nc ./INPUT_small/MOM.res.nc
ncks -O -d xaxis_1,0,71,2 -d yaxis_1,0,33,2 ./INPUT/cice.res.nc ./INPUT_small/cice.res.nc
