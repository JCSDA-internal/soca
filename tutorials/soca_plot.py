from netCDF4 import Dataset, Dataset
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
from datetime import datetime, timedelta
from dateutil import parser
import glob


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

# List the cycles
cycles=glob.glob('./2018*')
colors=['blue', 'red']
# Loop through obervation spaces and plot binned omb's
for k in obs_spaces:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for iter in [1, 2]:
        for c in cycles:
            cdate=datetime.strptime(c[-10:],'%Y%m%d%H')
            fnames=glob.glob(c+'/obs_out/outeriter_'+str(iter)+'/'+obs_spaces[k].ioname+'.out_00*.nc4')
            obs=Obs(fnames=fnames, varname=obs_spaces[k].varname)
            cycle_date='Initialized  \n'+cdate.strftime('%Y-%m-%d %H')
            obs.mae(ax, cdate, colors[iter-1])
    plt.title(obs_spaces[k].ioname, fontweight='bold')
    plt.xlabel('DA cycles', fontweight='bold')
    plt.ylabel('<|O-B|> '+obs_spaces[k].units, fontweight='bold')

    every_day = mdates.DayLocator(interval=1)
    ax.xaxis.set_major_locator(every_day)

    plt.grid(True)
    fig.savefig(k+'_iter_'+str(iter)+'.png')
