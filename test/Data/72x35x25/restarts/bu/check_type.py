import netCDF4 as nc

with nc.Dataset('MOM.res.nc', 'r') as ds:
    for dim_name, dim_obj in ds.dimensions.items():
        print(f"Dimension: {dim_name}, Size: {dim_obj.size}, {isinstance(dim_obj, float)}")
