#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Creates a configuration yaml file to write the analysis in a MOM6 restart
#

set -e

date_begin=$1
window_begin=$(date '+%Y-%m-%dT%H:%M:%S'Z -d "$date_begin")

cat << EOF > checkpoint.yaml

resolution:
  mom6_input_nml: ./inputnml/input.nml
  fields metadata: ./fields_metadata.yml

model:
  name: SOCA
  tstep: PT1H
  advance_mom6: 0
  model variables: [socn, tocn, hocn]
  tocn_minmax: [-1.8, 32.0]
  socn_minmax: [0.1, 38.0]

background:
  read_from_file: 1
  date: &date $window_begin
  basename: ./INPUT/
  ocn_filename: MOM.res.nc
  state variables: [socn, tocn, hocn]

analysis:
  read_from_file: 1
  date: *date
  basename: ./3dvar_out/
  ocn_filename: ocn.3dvarfgat.an.${window_begin}.nc
  state variables: [socn, tocn, hocn]
EOF
