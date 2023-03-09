#! /usr/bin/env python3

# (C) Copyright 2022-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
import shutil
from r2d2 import R2D2Data

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
workdir = conf['workdir']
if not os.path.exists(workdir):
    raise RuntimeError('Working directory does not exist')
os.chdir(workdir)

# Date
andate = conf['an']['date']
base = conf['experiment']['expid'] + '.an'  # + andate

# TODO put analysis back in restart file correctly?
file_types = ['MOM.res', 'cice.res']
expid = conf['experiment']['expid']
for in_pfx, out_pfx in (
        ('ice', 'cice.res'),
        ('ocn', 'MOM.res')):
    infile = f'{in_pfx}.{expid}.an.{andate}.nc'
    outfile = f'{expid}.an.{out_pfx}.nc'
    print(f"moving {infile} to {outfile}")
    shutil.move(infile, outfile)

for file_type in file_types:
    R2D2Data.store(
        model=conf['experiment']['model'],
        item='analysis',
        experiment=conf['experiment']['expid'],
        resolution=conf['resolution'],
        date=andate,
        domain='global',
        source_file=f'{base}.{file_type}.nc',
        file_extension='nc',
        file_type=file_type,
    )
