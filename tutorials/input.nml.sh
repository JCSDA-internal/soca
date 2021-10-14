#!/bin/bash
#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Template MOM6 namelist
#

set -e

year=$1
month=$2
day=$3
hour=$4

cat << EOF > input.nml
&MOM_input_nml
        output_directory = './',
        input_filename = 'r'
        restart_input_dir = 'INPUT/',
        restart_output_dir = 'RESTART/',
        parameter_filename = 'MOM_input' /

 &diag_manager_nml
 /

 &ocean_solo_nml
            months = 0
            days   = 1
            date_init = $year,$month,$day,$hour,0,0,
            hours = 0
            minutes = 0
            seconds = 0
            calendar = 'NOLEAP' /

 &fms_io_nml
      max_files_w=100
      checksum_required=.false.
/

 &fms_nml
       clock_grain='MODULE'
       domains_stack_size = 2000000
       clock_flags='SYNC' /
EOF
