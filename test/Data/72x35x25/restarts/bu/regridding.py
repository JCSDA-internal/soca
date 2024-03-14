import numpy as np
from scipy.interpolate import griddata
import netCDF4 as nc
from pyproj import Proj, transform

def regrid_latlon(source_lat, source_lon, target_lat, target_lon, data_values):
    # Convert lat-lon coordinates to Cartesian coordinates
    source_proj = Proj(proj='latlong', datum='WGS84')
    target_proj = Proj(proj='latlong', datum='WGS84')
    source_x, source_y = transform(source_proj, target_proj, source_lon, source_lat)

    target_x, target_y = np.meshgrid(target_lon, target_lat)
    
    # Interpolate data to the target grid
    data_regridded = griddata((source_x.flatten(), source_y.flatten()), data_values.flatten(),
                              (target_x, target_y), method='linear')

    return data_regridded

# Load the source NetCDF file with the original lat-lon coordinates
source_file_path = 'last_NOBM.res.nc'
source_dataset = nc.Dataset(source_file_path)

# Specify the target latitudes and longitudes
target_file_path = 'MOM.res.nc'
target_dataset = nc.Dataset(target_file_path)
target_lat = target_dataset.variables['lath'][:]  # Specify the target latitudes
target_lon = target_dataset.variables['lonh'][:]  # Specify the target longitudes

# Get the source latitudes, longitudes, and data values
source_lat = source_dataset.variables['lat'][:]
source_lon = source_dataset.variables['lon'][:]
data_values = source_dataset.variables['your_variable'][:]

# Perform regridding
regridded_data = regrid_latlon(source_lat, source_lon, target_lat, target_lon, data_values)

# Create a new NetCDF file for the regridded data
target_file_path = 'path/to/target_file.nc'
target_dataset = nc.Dataset(target_file_path, 'w', format='NETCDF4')

# Create dimensions in the target file
target_dataset.createDimension('lat', len(target_lat))
target_dataset.createDimension('lon', len(target_lon))

# Create variables in the target file
target_lat_var = target_dataset.createVariable('lat', 'f4', ('lat',))
target_lon_var = target_dataset.createVariable('lon', 'f4', ('lon',))
regridded_data_var = target_dataset.createVariable('your_variable_regridded', 'f4', ('lat', 'lon',))

# Assign data to variables
target_lat_var[:] = target_lat
target_lon_var[:] = target_lon
regridded_data_var[:, :] = regridded_data

# Close the datasets
source_dataset.close()
target_dataset.close()
