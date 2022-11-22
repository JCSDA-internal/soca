#!/bin/bash

# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set -eux

# cp ${SOCAREPO}/ewok/mom_input/MOM_input.${RESOLUTION}* ${WORKDIR}/MOM_input
cp -r ${SOCA_STATIC_DIR}/* ${WORKDIR}/
