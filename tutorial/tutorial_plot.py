from netCDF4 import Dataset, Dataset
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
from matplotlib import cm
from datetime import datetime, timedelta
from dateutil import parser
import glob
import sys
import cartopy.crs as ccrs
import cartopy

# State space classes and definitions
class Grid:
    def __init__(self,fname):
        ncfile = Dataset(fname,'r')
        self.lat=np.squeeze(ncfile.variables['lat'][:])
        self.lon=np.squeeze(ncfile.variables['lon'][:])
        self.mask=np.squeeze(ncfile.variables['mask2d'][:])
        ncfile.close()

def get_var(filename, varname):
    ncfile = Dataset(filename,'r')
    var = np.squeeze(ncfile.variables[varname][:])
    ncfile.close()
    return var

def plothor(x, y, z, varname='', clim=[-1,1], pngname='', title=''):
    fig = plt.figure(figsize=(13,6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    clevs = np.linspace(clim[0], clim[1], 41)
    p = ax.contourf(x, y, z, clevs,
                    transform=ccrs.PlateCarree(),
                    extend='both',
                    cmap=cm.bwr)

    plt.colorbar(p, ax=ax, shrink=0.5)
    plt.title(title)
    ax.coastlines()
    fig.savefig(pngname)

# Observation space classes and definitions
class Obs():
    def __init__(self, fnames="",varname=""):
        self.fnames = fnames
        self.varname = varname
        self.lon = []
        self.lat = []
        self.time = []
        self.omb = []
        self.obs = []
        self.hofx = []
        self.qc = []

        for f in fnames:
            self._read(f)

        for t in range(len(self.time)):
            self.time[t]=parser.parse(self.time[t]).timestamp()/3600.0

    def _read(self, fname):
        ncfile = Dataset(fname,'r')
        grp = ncfile.groups['MetaData']
        self.lon = np.append(self.lon, grp.variables['longitude'][:])
        self.lat = np.append(self.lat, grp.variables['latitude'][:])
        self.time = np.append(self.time, grp.variables['datetime'][:])

        grp = ncfile.groups['ombg']
        self.omb = np.append(self.omb, -grp.variables[self.varname][:])

        grp = ncfile.groups['ObsValue']
        self.obs = np.append(self.obs, grp.variables[self.varname][:])

        grp = ncfile.groups['hofx']
        self.hofx = np.append(self.hofx, grp.variables[self.varname][:])

        grp = ncfile.groups['EffectiveQC1']
        self.qc = np.append(self.qc, grp.variables[self.varname][:])

    def mae(self, ax, cycle_date, color):
        qci = np.where(self.qc==0)
        time = self.time[qci]
        omb = np.abs(self.omb[qci])
        mae = []
        time_range, dt = np.linspace(np.min(time),np.max(time),24, retstep=True)
        for t in time_range:
            i=np.where( np.logical_and(time>t-0.5*dt, time<=t+0.5*dt) )
            mae.append(np.mean(np.abs(omb[i])))
        td = []
        for t in range(24):
            td.append( cycle_date + timedelta(hours=t) )
        ax.plot_date(td, mae, '-', color=color)

# Define the observation space dictionary
class ObsSpace():
    def __init__(self, ioname, varname, units):
        self.ioname = ioname
        self.varname = varname
        self.units = units

obs_spaces={
    "sst": ObsSpace("sst", "sea_surface_temperature","[K]"),
    "sss": ObsSpace("sss", "sea_surface_salinity", "[psu]"),
    "adt": ObsSpace("adt", "obs_absolute_dynamic_topography", "[m]"),
    "T": ObsSpace("insitu.T", "sea_water_temperature", "[K]"),
    "S": ObsSpace("insitu.S", "sea_water_salinity", "[psu]")
}


if __name__ == '__main__':

    # get the experiment type from the comman line
    plot_type = sys.argv[1]

    # read the grid
    grid = Grid(fname='../static/soca_gridspec.nc')

    if (plot_type == '3dvar'):
        fname = './3dvar_out/ocn.3dvar.incr.2018-04-15T00:00:00Z.nc'

        # ssh increment
        ssh = get_var(filename = fname, varname = 'ave_ssh')
        plothor(grid.lon, grid.lat, ssh, varname='ssh',
                clim=[-0.2,0.2], pngname='incr.ssh.png', title='SSH increment [m]')

        # sst increment
        sst = get_var(filename = fname, varname = 'Temp')
        plothor(grid.lon, grid.lat, sst[0,:,:],
                varname='sst', clim=[-5,5], pngname='incr.sst.png', title='SST increment [k]')

        # sst increment
        sss = get_var(filename = fname, varname = 'Salt')
        plothor(grid.lon, grid.lat, sss[0,:,:],
                varname='sss', clim=[-1,1], pngname='incr.sss.png', title='SSS increment [psu]')

    if (plot_type == '3dvarfgat'):
        # List the cycles
        cycles=glob.glob('./2018*')
        colors=['blue', 'red']

        # Loop through obervation spaces and plot time series of global mean observation - background
        for k in obs_spaces:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for iter in [1, 2]:
                for c in cycles:
                    cdate = datetime.strptime(c[-10:],'%Y%m%d%H')
                    fnames = glob.glob(c+'/obs_out/outeriter_'+str(iter)+'/'+obs_spaces[k].ioname+'.out_00*.nc4')
                    obs = Obs(fnames=fnames, varname=obs_spaces[k].varname)
                    cycle_date = 'Initialized  \n'+cdate.strftime('%Y-%m-%d %H')
                    obs.mae(ax, cdate, colors[iter-1])

            ax.legend([ax.lines[0], ax.lines[len(cycles)]],['outer iteration 1', 'outer iteration 2'])
            plt.title(obs_spaces[k].ioname, fontweight='bold')
            plt.xlabel('DA cycles', fontweight='bold')
            plt.ylabel('<|O-B|> '+obs_spaces[k].units, fontweight='bold')
            every_day = mdates.DayLocator(interval=1)
            ax.xaxis.set_major_locator(every_day)
            plt.grid(True)
            fig.savefig(k+'_global_mae.png')

        # Outer loop increments for the first cycle
        for iter in [1, 2]:
            fname = './2018041700/ocn.3dvarfgat.iter'+str(iter)+'.incr.2018-04-17T00:00:00Z.nc'

            # ssh increment
            ssh = get_var(filename = fname, varname = 'ave_ssh')
            plothor(grid.lon, grid.lat, ssh,
                    varname='ssh', clim=[-0.2,0.2], pngname='incr.'+str(iter)+'.ssh.png',
                    title='outer iteration '+str(iter)+' SSH increment [m]')

            # sst increment
            sst = get_var(filename = fname, varname = 'Temp')
            plothor(grid.lon, grid.lat, sst[0,:,:],
                    varname='sst', clim=[-5,5], pngname='incr.'+str(iter)+'.sst.png',
                    title='outer iteration '+str(iter)+' SST increment [k]')

            # sss increment
            sss = get_var(filename = fname, varname = 'Salt')
            plothor(grid.lon, grid.lat, sss[0,:,:],
                    varname='sss', clim=[-1,1], pngname='incr.'+str(iter)+'.sss.png',
                    title='outer iteration '+str(iter)+' SSS increment [psu]')
