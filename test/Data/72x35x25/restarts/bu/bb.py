import numpy as np
#from scipy.interpolate import griddata
import netCDF4 as nc
from pyproj import Proj, transform


import xarray as xr

# Replace 'input_file.nc' with the path to your NetCDF file
ds = xr.open_dataset('bu_NOBM.res.nc')

target_file_path = 'MOM.res.nc'
target_dataset = nc.Dataset(target_file_path)
target_lats = target_dataset.variables['lath'][:]  # Specify the target latitudes                                                       
target_lons = target_dataset.variables['lonh'][:]  # Specify the target longitudes  

# Replace 'target_lats' and 'target_lons' with your tripolar grid coordinates
#target_lats = [...]  # List or array of target latitudes
#target_lons = [...]  # List or array of target longitudes

target_grid = xr.DataArray(0, coords=[('lat', target_lats), ('lon', target_lons)])

# Replace 'lat' and 'lon' with the original latitude and longitude coordinates in your dataset
ds_regridded = ds.interp(lat=target_lats, lon=target_lons, method='linear')

# Replace 'output_file.nc' with the desired output file name
ds_regridded.to_netcdf('output_file.nc')
