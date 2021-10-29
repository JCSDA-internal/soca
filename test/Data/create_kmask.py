#!/usr/bin/env python3
import netCDF4
import numpy as np

# Invent the Ice/Ocean Jacobian: dc/dt
ncf=netCDF4.Dataset("./INPUT/cice.res.nc",'r')
aice_mask=ncf.variables['aicen'][:]
ncf.close()
dcdt=-0.01*np.ones((35,72))*aice_mask

# Convert the dirac output into a mask for the dynamic height balance
ncf=netCDF4.Dataset("./Data/ocn.balance_mask.an.2018-04-15T00:00:00Z.nc",'r')
kmask=ncf.variables['ave_ssh'][:]
ncf.close()
kmask[kmask<0.0]=0.0
kmask=1.0-kmask

# Write to file
ncf=netCDF4.Dataset("./Data/kmask.nc",'w')
nx=ncf.createDimension('xaxis_1',72)
ny=ncf.createDimension('yaxis_1',35)
v1=ncf.createVariable('dcdt','f8',('yaxis_1','xaxis_1'))
v2=ncf.createVariable('kmask','f8',('yaxis_1','xaxis_1'))
v1[:]=dcdt
v2[:]=kmask
ncf.close()
