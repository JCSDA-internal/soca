#!/usr/bin/env python
import numpy as np
import netCDF4 as nc
from scipy.spatial import KDTree


vars = ['liq_precip','froz_precip','liq_runoff', 'froz_runoff']
s_times = [3]
src_file = '360x210x75/INPUT/ocean_precip_monthly.nc'
dst_file = '90x43x25/INPUT/forcing_monthly.nc'

# vars = ['SW', 'LW', 'latent', 'sensible', 'evap', 'taux','tauy', 'SST', 'SSS']
# s_times = [103,104,105]
# src_file = '360x210x75/INPUT/ocean_forcing_daily.nc'
# dst_file = '90x43x25/INPUT/forcing_daily.nc'

# load source grid
d='/home/tsluka/work/jedi-soca/mom-configs/360x210x75.hybrid'
ncd=nc.Dataset(d+'/ocean_geometry.nc','r')
src_lat = ncd.variables['geolat'][:]
src_lon = ncd.variables['geolon'][:]
src_lath = ncd.variables['lath'][:]
src_lonh = ncd.variables['lonh'][:]
src_msk = ncd.variables['wet'][:] > 0.5


# load destination grid
d='/home/tsluka/work/jedi-soca/mom-configs/90x43x25.hybrid'
ncd=nc.Dataset(d+'/ocean_geometry.nc','r')
dst_lat = ncd.variables['geolat'][:]
dst_lon = ncd.variables['geolon'][:]
dst_msk = ncd.variables['wet'][:] > 0.5

# create kd-tree lookup for source grid
src_x, src_y = np.meshgrid(np.arange(src_msk.shape[1]), np.arange(src_msk.shape[0]))
src_x = src_x[src_msk].ravel()
src_y = src_y[src_msk].ravel()
src_lat = src_lat[src_msk].ravel()
src_lon = src_lon[src_msk].ravel()
kd = KDTree(zip(src_lon, src_lat))

# open source file
s_ncd = nc.Dataset(src_file)

# create output file
d_ncd = nc.Dataset(dst_file, 'w')
d_ncd.createDimension('xh', dst_msk.shape[1])
d_ncd.createDimension('yh', dst_msk.shape[0])
d_ncd.createDimension('time')
for v in vars:
    d_ncd.createVariable(v, 'f4', ('time','yh','xh'))
d_ncd.createVariable('time', 'f8', ('time',))
d_ncd.createVariable('xh', 'f8', ('xh',))
d_ncd.createVariable('yh', 'f8', ('yh',))

# create output data
for v in vars:
    val_d = np.zeros( (len(s_times), dst_msk.shape[0], dst_msk.shape[1]) )
    val_s = s_ncd.variables[v][:]
    for y in range(dst_msk.shape[0]):
        for x in range(dst_msk.shape[1]):
            if dst_msk[y,x]:
                _, n = kd.query( (dst_lon[y,x], dst_lat[y,x]))
                for td, ts in enumerate(s_times):
                    val_d[td,y,x] = val_s[ts, src_y[n], src_x[n] ]
            else:
                val_d[:,y,x] = 0.0
    d_ncd.variables[v][:] = val_d

# fill in time
val_d = d_ncd.variables['time'][:]
val_s = s_ncd.variables['time'][:]
for td, ts in enumerate(s_times):
    val_d[td] = val_s[ts]
d_ncd.variables['time'][:] = val_d
for k in s_ncd.variables['time'].ncattrs():
    d_ncd.variables['time'].setncattr( k, s_ncd.variables['time'].getncattr(k))
d_ncd.variables['xh'][:] = ncd.variables['lonh'][:]
d_ncd.variables['yh'][:] = ncd.variables['lath'][:]
for k in ncd.variables['lonh'].ncattrs():
    d_ncd.variables['xh'].setncattr( k, ncd.variables['lonh'].getncattr(k))
for k in ncd.variables['lath'].ncattrs():
    d_ncd.variables['yh'].setncattr( k, ncd.variables['lath'].getncattr(k))


# done
d_ncd.grid_type="regular"
d_ncd.close()

