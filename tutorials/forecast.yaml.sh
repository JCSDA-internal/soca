#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Creates a 24 hour forecast configuration yaml file
#

set -e

date_begin=$1
window_begin=$(date '+%Y-%m-%dT%H:%M:%S'Z -d "$date_begin")

cat << EOF > forecast.yaml

forecast length: P1D

geometry:
  mom6_input_nml: ./inputnml/input.nml
  fields metadata: ./fields_metadata.yml

model:
  name: SOCA
  tstep: PT1H
  advance_mom6: 1
  model variables: &soca_vars [socn, tocn, ssh, hocn]

initial condition:
  read_from_file: 1
  date: &date $window_begin
  basename: ./INPUT/
  ocn_filename: MOM.res.nc
  state variables:  *soca_vars

output:
  frequency: PT1H
  datadir: fcst
  exp: mom6
  date: *date
  type: fc
EOF
