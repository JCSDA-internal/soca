import xarray as xr
import xesmf as xe
import netCDF4 as nc

# Load the source NetCDF file
source_file_path = 'bu_NOBM.res.nc'
source_ds = xr.open_dataset(source_file_path)

# Create a grid for the target lat-lon coordinates

target_file_path = 'MOM.res.nc'
target_dataset = nc.Dataset(target_file_path)
target_lats = target_dataset.variables['lath'][:]  # Specify the target latitudes                                                                             
target_lons = target_dataset.variables['lonh'][:]  # Specify the target longitudes  
target_grid = {'lat': target_lats, 'lon': target_lons}

# Create regridder
regridder = xe.Regridder(source_ds, target_grid, method='bilinear', periodic=True)

# Perform the regridding
target_ds = regridder(source_ds)

# Save the result to a new NetCDF file
target_file_path = 'remapped_file.nc'
target_ds.to_netcdf(target_file_path)
