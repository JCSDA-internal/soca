#!/usr/bin/env python

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as stats
import cartopy.crs as ccrs
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from cmocean import cm as cmo

def plothor(x, y, z, map, varname='', clim=[-1,1], proj_type='reg', plot_type='pcolor', colormap=cmo.balance):
    a=clim[0]
    b=clim[1]

    if  proj_type == 'reg':
        proj = ccrs.Mollweide()
    if proj_type == 'north':
        proj = ccrs.NorthPolarStereo()

    fig = plt.figure(figsize=(13,6.2))
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    if  proj_type == 'reg':
        ax.set_global()
    if proj_type == 'north':
        ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

    if plot_type == 'contour':
        clevs = np.linspace(a, b, 41)
        p = ax.contourf(x, y, z, clevs,
                        transform=ccrs.PlateCarree(),
                        cmap=colormap,
                        extend='both')
    if plot_type == 'pcolor':
        p = ax.pcolormesh(x, y, z,
                          vmin=clim[0],
                          vmax=clim[1],
                          transform=ccrs.PlateCarree(),
                          cmap=colormap )

    plt.colorbar(p, ax=ax)
    plt.title(varname)
    ax.coastlines()

class Grid:
    def __init__(self,fname):
        ncfile = Dataset(fname,'r')
        self.lat=np.squeeze(ncfile.variables['lat'][:])
        self.lon=np.squeeze(ncfile.variables['lon'][:])
        self.mask=np.squeeze(ncfile.variables['mask2d'][:])
        ncfile.close()

def get_var(filename, varname, level=0):
    ncfile = Dataset(filename,'r')
    tmp=np.squeeze(ncfile.variables[varname][:])
    ncfile.close()
    try:
        tmp=tmp[level,:,:]
    except:
        pass
    return tmp

if __name__ == '__main__':
    description = """Plot horizontal field, global and north polar stereo:
                     plt_horiz fld.py -g /path/to/socagrid/soca_gridspec.nc
                                      -f /path/to/socafile/incr.1.nc
                                      -v ssh
                                      -l -1 1
                  """

    parser = ArgumentParser(
        description=description,
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-g',
        '--grid',
        help="soca geometry file, soca_gridspec.nc",
        type=str, required=True)
    parser.add_argument(
        '-f',
        '--file',
        help="soca file containg fields",
        type=str, required=True)
    parser.add_argument(
        '-v',
        '--variable',
        help="Variable name",
        type=str, required=True)
    parser.add_argument(
        '-l',
        '--clim',
        help="Color scale: [min, max]",
        type=float, nargs='+', required=True)
    parser.add_argument(
        '-c',
        '--color',
        help="Color map name from cmocean: cmo.balance, cmo.thermal, cmo.haline, ...",
        type=str, required=True)

    args = parser.parse_args()
    grid = Grid(fname=args.grid)
    var = get_var(filename=args.file, varname=args.variable)
    plothor(grid.lon,grid.lat,var/grid.mask,
            map,varname=args.variable,
            clim=args.clim,
            colormap=args.color)
    plothor(grid.lon,grid.lat,var/grid.mask,map,
            varname=args.variable,
            clim=args.clim,
            proj_type='north',
            colormap=args.color)
    plt.show()
