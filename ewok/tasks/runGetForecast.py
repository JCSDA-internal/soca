#! /usr/bin/env python3

# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
from r2d2 import R2D2Data

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['workdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['workdir'])

# Define experiment to read from, current experiment by default
exp_read = conf['experiment']['expid']
if 'exp_source' in conf:
    exp_read = conf['exp_source']

# Fetch state
base = conf['experiment']['expid'] + '.fc.'
sdate = conf['fcdate'] + '.' + conf['fcstep']
# filename = base + sdate + '.$(file_type).nc'

# determine the files needed depending on resolution, and if ice is on or off
file_types = ['MOM.res', 'cice.res']
if conf['resolution'] in {'025deg',}:
    file_types += ['MOM.res_1', 'MOM.res_2', 'MOM.res_3']
elif conf['resolution'] in {'1deg','5deg'}:
    pass
else:
    raise RuntimeError(f'Unsupported resolution {conf["resolution"]}')

for file_type in file_types:
    R2D2Data.fetch(
        model=conf['experiment']['model'],
        item='forecast',
        experiment=exp_read,
        step=conf['fcstep'],
        resolution=conf['resolution'],
        date=conf['fcdate'],
        domain='global',
        target_file=f'{base}{sdate}.{file_type}.nc',
        file_extension='nc',
        file_type=file_type,
    )
