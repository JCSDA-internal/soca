#
# (C) Copyright 2021-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Create ioda file with random locations.
#

from netCDF4 import Dataset, num2date, date2num
import numpy as np
import datetime

class Ioda():
    def __init__(self, varname, fname_out, datetime_start, datetime_end, nobs, sigo, dodepth=False):
        self.datetime = []
        self.lon = []
        self.lat = []
        self.obserr = sigo
        self.depth = []
        self.dodepth = dodepth
        self.nobs = nobs
        self.varname=varname
        self.time_start=datetime_start
        self.time_end=datetime_end
        self.fname_out=fname_out

        # Create random locations
        self._synth_obs()

        # Save a minimalist ioda file
        self._obs2ioda()

    def _synth_obs(self):
        Dt=self.time_end-self.time_start
        dts = np.random.uniform(low=0.0, high=Dt.total_seconds(), size=(self.nobs,))

        for dt in dts:
            timeiso=(self.time_start+datetime.timedelta(seconds=int(dt))).isoformat()+'Z'
            self.datetime=np.append(self.datetime, timeiso)

        self.datetime=np.array(self.datetime)
        self.lon=np.random.uniform(low=-180.0, high=180.0, size=(self.nobs,))
        self.lat=np.random.uniform(low=-90.0, high=90.0, size=(self.nobs,))
        if self.dodepth:
            self.depth=np.random.uniform(low=0.0, high=2000.0, size=(self.nobs,))

    def _obs2ioda(self):
        ncfile = Dataset(self.fname_out,'w',format='NETCDF4')

        # Dimensions
        nlocs_dim = ncfile.createDimension('nlocs', None)

        # Global attributes
        ncfile.nvars=1
        ncfile.nlocs = len(self.datetime)

        # Variables
        tmp = ncfile.createVariable('nlocs', np.int32, ('nlocs',))
        tmp[:] = np.arange(1, len(self.datetime)+1)

        # Group MetaData
        grp = ncfile.createGroup('MetaData')
        tmp = grp.createVariable('datetime',str,('nlocs',))
        tmp[:] = self.datetime
        tmp = grp.createVariable('latitude',np.float32,('nlocs',))
        tmp[:] = self.lat
        tmp = grp.createVariable('longitude',np.float32,('nlocs',))
        tmp[:] = self.lon
        if self.dodepth:
            tmp = grp.createVariable('depth',np.float32,('nlocs',))
            tmp[:] = self.depth

        # Group ObsError
        grp = ncfile.createGroup('ObsError')
        tmp = grp.createVariable(self.varname,np.float32,('nlocs',))
        tmp[:]=self.obserr*np.ones(np.shape(self.datetime))

        # Group ObsValue
        grp = ncfile.createGroup('ObsValue')
        tmp = grp.createVariable(self.varname,np.float32,('nlocs',))
        tmp[:] = -9999.9*np.ones(np.shape(self.datetime))

        # Group PreQC
        grp = ncfile.createGroup('PreQC')
        tmp = grp.createVariable(self.varname,np.int32,('nlocs',))
        preqc = 0*np.zeros(np.shape(self.datetime))
        tmp[:] = preqc

        ncfile.close()

time_start = datetime.datetime(2018, 4, 15)
time_end = datetime.datetime(2018, 4, 20)

# Sea surface temperature
Ioda(varname='sea_surface_temperature',
     fname_out='./obs_scratch/sst.nc4',
     datetime_start=time_start,
     datetime_end=time_end,
     sigo=0.5,
     nobs=10000)

# Sea surface salinity
Ioda(varname='sea_surface_salinity',
     fname_out='./obs_scratch/sss.nc4',
     datetime_start=time_start,
     datetime_end=time_end,
     sigo=0.2,
     nobs=10000)

# Absolute dynamic topography
Ioda(varname='obs_absolute_dynamic_topography',
     fname_out='./obs_scratch/adt.nc4',
     datetime_start=time_start,
     datetime_end=time_end,
     sigo=0.1,
     nobs=8000)

# Insitu temperature
Ioda(varname='sea_water_temperature',
     fname_out='./obs_scratch/insitu.T.nc4',
     datetime_start=time_start,
     datetime_end=time_end,
     sigo=0.5,
     nobs=1000, dodepth=True)

# Insitu salinity
Ioda(varname='sea_water_salinity',
     fname_out='./obs_scratch/insitu.S.nc4',
     datetime_start=time_start,
     datetime_end=time_end,
     sigo=0.2,
     nobs=500, dodepth=True)
