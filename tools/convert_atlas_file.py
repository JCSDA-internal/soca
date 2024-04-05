#!/usr/bin/env python3
"""
Tool to convert from the ATLAS unstructured output, to the same grid used by SOCA/MOM6.
Useful for visualizing output from ATLAS
"""
import netCDF4 as nc
import numpy as np
from scipy.spatial import KDTree
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('inputFilename', type=str, help='The input file name')
parser.add_argument('outputFilename', type=str, help='The output file name')
parser.add_argument('--grid', type=str, default='./soca_gridspec.nc',
                    help='The output grid file name (default: ./soca_gridspec.nc)')
args = parser.parse_args()

outputGrid = args.grid
inputFileName = args.inputFilename
outputFileName = args.outputFilename

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

  # find variables to process
  allVars = [name for name, var in ncdIn.variables.items() if len(var.dimensions) == 2]
  for v in allVars:
    valSrc = ncdIn.variables[v]
    nz = valSrc.shape[-1]
    zName = f'nz_{v}'
    ncdOut.createDimension(zName, nz)
    valDst = ncdOut.createVariable(v, valSrc.datatype, [zName,'latitude','longitude'])
    valDst[:] = valSrc[:].reshape((ny,nx,nz)).transpose(2,0,1)
