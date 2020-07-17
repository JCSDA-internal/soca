#!/usr/bin/env python
import netCDF4
import numpy as np

dcdt=-0.01*np.ones((35,72))

ncf=netCDF4.Dataset("./Data/dcdt.nc",'w')
nx=ncf.createDimension('xaxis_1',72)
ny=ncf.createDimension('yaxis_1',35)
v=ncf.createVariable('dcdt','f8',('yaxis_1','xaxis_1'))
v[:]=dcdt
ncf.close()
