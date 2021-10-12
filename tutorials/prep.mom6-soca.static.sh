#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Copy or link the MOM6 and soca static files
#

set -e

datadir=$1

# Prepare soca/MOM6
mkdir -p RESTART
cp -rL $datadir/Data/72x35x25/INPUT .
mkdir -p ./inputnml
cp $datadir/Data/72x35x25/input.nml ./inputnml/
cp $datadir/Data/72x35x25/MOM_input .
ln -sf $datadir/Data/fields_metadata.yml .
ln -sf $datadir/Data/72x35x25/*_table .
ln -sf $datadir/Data/rossrad.dat .
ln -sf $datadir/Data/godas_sst_bgerr.nc .
