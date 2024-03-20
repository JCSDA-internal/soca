#!/usr/bin/env python3
"""
Tool to convert from the ATLAS unstructured output, to the same grid used by SOCA/MOM6.
Useful for visualizing output from ATLAS
"""
import netCDF4 as nc
import numpy as np
from scipy.spatial import KDTree

outputGrid = "../gridgen/soca_gridspec.72x35x25.nc"
inputFileName = "parametricStdDev.nc"
outputFileName = "output.nc"

# get destination grid coordinates
with nc.Dataset(outputGrid, 'r') as ncd:
  dstLon = ncd.variables['lon'][0,:,:]
  dstLat = ncd.variables['lat'][0,:,:]
  datMask = ncd.variables['mask2d'][0,:,:]
  ny, nx = dstLon.shape
dstPoints = np.column_stack((dstLat.flatten(), dstLon.flatten()))

with nc.Dataset(inputFileName, 'r') as ncdIn, nc.Dataset(outputFileName, 'w', format='NETCDF4') as ncdOut:
  # initialize output file lat/lon
  ncdOut.createDimension('latitude', ny)
  ncdOut.createDimension('longitude', nx)
  latitudes = ncdOut.createVariable('latitude', 'f4', ('latitude',))
  longitudes = ncdOut.createVariable('longitude', 'f4', ('longitude',))  
  latitudes[:] = np.mean(dstLat, axis=1)  # Assuming dstLat is 2D
  longitudes[:] = np.mean(dstLon, axis=0)  # Assuming dstLon is 2D  

  # Assign attributes to latitude and longitude variables
  latitudes.units = 'degrees_north'
  latitudes.standard_name = 'latitude'
  longitudes.units = 'degrees_east'
  longitudes.standard_name = 'longitude'

  # # use KD tree to find mapping between src and dst points
  # srcLon = ncdIn.variables['lon'][:]
  # srcLat = ncdIn.variables['lat'][:]
  # srcPoints = np.column_stack((srcLat.flatten(), srcLon.flatten()))
  # kdTree = KDTree(srcPoints)
  # _, mapping = kdTree.query(dstPoints)

  # find variables to process
  allVars = [name for name, var in ncdIn.variables.items() if len(var.dimensions) == 2]
  for v in allVars:
    valSrc = ncdIn.variables[v]
    nz = valSrc.shape[-1]
    zName = f'nz_{v}'
    ncdOut.createDimension(zName, nz)
    valDst = ncdOut.createVariable(v, valSrc.datatype, [zName,'latitude','longitude'])
    #valDst[:] = valSrc[mapping, :].reshape((ny,nx,nz)).transpose(2,0,1)
    valDst[:] = valSrc[:].reshape((ny,nx,nz)).transpose(2,0,1)
