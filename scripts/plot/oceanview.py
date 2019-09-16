#!/usr/bin/env python

from __future__ import print_function
import matplotlib as mpl
mpl.use('WXAgg')
mpl.interactive(False)
import pylab as pl
from pylab import get_current_fig_manager as gcfm
import wx
import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
from tqdm import tqdm
from scipy import stats
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def get_from_ioda(ncfile, varname):
    varout=ncfile.variables[varname][:]
    return varout

def plot_prof(dax, fcst, ana, obs, z, var_name):
    dax.plot(fcst, z, '-og',lw=1)
    dax.plot(ana, z, '-sr',lw=2)
    dax.plot(obs, z, '-*b')
    dax.set_xlabel(var_name)
    dax.grid(True)

def plot_reg(dax, x1, x2, xcolor, title):
    dax.plot(x1, x2,'.', color=xcolor, alpha=0.1)
    dax.axis('equal')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x1, x2)
    titlestr=title + " slope: %3.3f r: %1.3f rms: %1.3f" % (slope, r_value, std_err)
    plt.grid(True)
    plt.title(titlestr,fontsize=24,fontweight='bold')

def draw_map(lonl=-180, lonr=180):
    map = Basemap(projection='robin',lon_0=0,resolution='c')
    map.drawcoastlines()
    map.fillcontinents(color='brown',lake_color='aqua', zorder=1)
    map.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
    map.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1])
    map.drawmapboundary(fill_color='aqua',zorder=0)
    return map

def find_inst(INSTID, dict_inst):
    for instrument in dict_inst:
        if INSTID==dict_inst[instrument].instid:
            inst_name=dict_inst[instrument].name
    return inst_name

def find_var(VARID, dict_varspecs):
    for var in dict_varspecs:
        if VARID==dict_varspecs[var].varid:
            var_name=var
    return var_name

def obsid2instidvarid(obsid, var):
    instid=np.ones(np.shape(obsid))
    varid=np.ones(np.shape(obsid))
    if (var=='sea_water_temperature'):
        varid=101*varid
        instid=508*instid
    if (var=='sea_water_salinity'):
        varid=102*varid
        instid=508*instid
    if (var=='obs_absolute_dynamic_topography'):
        varid=105*varid
        instid=517*instid
    if (var=='sea_surface_temperature'):
        varid=101*varid
        instid=520*instid
    return instid, varid

class Instrument:
    def __init__(self, name, instid, varid, zmin, zmax, color='k'):
        self.name = name
        self.instid = instid      # Instrument ID
        self.varid = varid        # Var ID
        self.zmin = zmin
        self.zmax = zmax
        self.color = color

class VarSpecs:
    def __init__(self, varid, bounds, fignum, units):
        self.varid  = varid
        self.bounds = bounds
        self.fignum = fignum
        self.units  = units

def dictionaries():
    dict_inst = {
                 'PROFILE'   : Instrument(name='PROFILE',   instid=508, varid=np.array([101, 102]), zmin=0, zmax=2000),
                 'ALTIMETER' : Instrument(name='ALTIMETER', instid=517, varid=np.array([105]),      zmin=0, zmax=0),
                 'SST'       : Instrument(name='SST',       instid=520, varid=np.array([101]),      zmin=0, zmax=10),
    }

    dict_varspecs = {'T'    : VarSpecs(varid=101, bounds= [-3., 38.], fignum=1, units='^oC'),
                     'S'    : VarSpecs(varid=102, bounds= [ 0., 45.], fignum=2, units='psu'),
                     'U'    : VarSpecs(varid=103, bounds= [-5.,  5.], fignum=3, units='m/s'),
                     'V'    : VarSpecs(varid=104, bounds= [-5.,  5.], fignum=4, units='m/s'),
                     'SSH'  : VarSpecs(varid=105, bounds= [-4.,  4.], fignum=5, units='m'),
                     'aice' : VarSpecs(varid=106, bounds= [0.,   1.], fignum=6, units=''),
                     'hice' : VarSpecs(varid=107, bounds= [0.,  20.], fignum=7, units='m')}

    COLORS=['b','r','g','k','m','y']

    return dict_inst, dict_varspecs, COLORS

class ioda:
    def __init__(self, iodafnames):
        flist=iodafnames
        self.lon=[]
        self.lat=[]
        self.oma=[]
        self.omf=[]
        self.obs=[]
        self.obserror=[]
        self.lev=[]
        self.varid=[]
        self.instid=[]
        self.col=[]
        self.tt=[]
        list_of_vars = ['sea_water_temperature',
                        'sea_water_salinity',
                        'obs_absolute_dynamic_topography',
                        'sea_surface_temperature']
        for iodafname in tqdm(flist):
            ncfile = Dataset(iodafname)
            for var in list_of_vars:
                try:
                    try:
                        ivar=var+'@ombg'
                        dum = get_from_ioda(ncfile,ivar)
                        I=np.where(abs(dum)<9999.9)
                        self.omf=np.append(-dum[I],self.omf)
                    except:
                        ivar=var+'@ombg1'
                        dum = get_from_ioda(ncfile,ivar)
                        I=np.where(abs(dum)<9999.9)
                        self.omf=np.append(-dum[I],self.omf)

                    ivar=var+'@oman'
                    dum=get_from_ioda(ncfile,ivar);
                    self.oma = np.append(-dum[I],self.oma)
                    ivar=var+'@ObsValue'
                    dum=get_from_ioda(ncfile,ivar);
                    self.obs = np.append(dum[I],self.obs)
                    ivar=var+'@ObsError'
                    dum=get_from_ioda(ncfile,ivar);
                    self.obserror = np.append(dum[I],self.obserror)

                    dum=get_from_ioda(ncfile,'longitude@MetaData');
                    self.lon = np.append(dum[I],self.lon)
                    dum=get_from_ioda(ncfile,'latitude@MetaData');
                    self.lat = np.append(dum[I],self.lat)
                    try:
                        dum=get_from_ioda(ncfile,'depth@MetaData');
                        self.lev = np.append(-dum[I],self.lev)
                    except:
                        self.lev = np.append(0*dum[I],self.lev)
                    instid, varid = obsid2instidvarid(dum, var)
                    self.instid = np.append(instid, self.instid)
                    self.varid = np.append(varid, self.varid)
                except:
                    pass


            ncfile.close()
        nobs=len(self.lon)
        self.tt = np.zeros(np.shape(self.varid))

class observation_space(object):
    def __init__(self,iodafname,i=[]):
        self.ioda=ioda(iodafname)
        self.iodafname=iodafname
        self.fignum = 2
        self.figure = pl.figure(num=1, figsize=(18, 10))
        self.axis = self.figure.add_subplot(111)
        self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100))
        gcfm().canvas.SetToolTip(self.tooltip)
        self.tooltip.Enable(False)
        self.tooltip.SetDelay(0)
        self.figure.canvas.mpl_connect('motion_notify_event', self._onMotion)
        self.figure.canvas.mpl_connect('button_press_event', self._onClick)
        self.dataX = np.squeeze(self.ioda.lon)
        self.dataY = np.squeeze(self.ioda.lat)

        map0=draw_map(lonl=-180,lonr=180)
        x, y =map0(self.dataX,self.dataY)
        self.X=x
        self.Y=y
        cnt=0
        for inst in [508, 517, 520]:
            I=np.where(self.ioda.instid==inst)
            self.axis.plot(x[I], y[I], linestyle='None', marker='.', markersize=5, label='myplot',color=COLORS[cnt])
            cnt+=1

    def _onMotion(self, event):
        collisionFound = False
        if event.xdata != None and event.ydata != None: # mouse is inside the axes
            I=np.where( (self.ioda.instid!=520)|(self.ioda.instid!=517) )
            for i in range(len(self.X)):
                radius = 100000 # Collision radius
                if (abs(event.xdata - self.X[i]) < radius) and (abs(event.ydata - self.Y[i]) < radius):
                    inst_name = find_inst(self.ioda.instid[i], dict_inst)
                    var_name = find_var(self.ioda.varid[i], dict_varspecs)
                    top = tip='Lon=%f\nLat=%f\nInstrument: %s\nVar: %s' % (self.dataX[i], self.dataY[i], inst_name, var_name)
                    self.tooltip.SetTip(tip)
                    self.tooltip.Enable(True)
                    self.i=i
                    collisionFound = True
                    break
        if not collisionFound:
            self.tooltip.Enable(False)

    def _onClick(self, event):

        if event.button == 1:    # Left mouse click: Profile
            self.figure2 = plt.figure(num=self.fignum, figsize=(16, 12), facecolor='c')
            self.axis2 = self.figure2.add_axes([0.3,0.69,0.4,0.3])
            map=draw_map()
            for shift in [0, 360]:
                x, y =map(self.dataX+shift,self.dataY)
                self.axis2.plot(x[:], y[:], linestyle='None', marker='.', markersize=2, label='myplot',color='b')
                self.axis2.plot(x[self.i], y[self.i], linestyle='None', marker='.', markersize=10, label='myplot', color='k')

            # Identify instrument and variable
            inst_name = find_inst(self.ioda.instid[self.i], dict_inst)
            var_name = find_var(self.ioda.varid[self.i], dict_varspecs)

            # Prepare axis
            self.axis3 = self.figure2.add_axes([0.1,0.05,0.35,0.6])
            self.axis3.set_ylabel('Depth [m]')
            self.axis4 = self.figure2.add_axes([0.55,0.05,0.35,0.6])

            # Get indices of observation pointed by mouth
            I=np.where( (self.ioda.lon==self.dataX[self.i]) & (self.ioda.lat==self.dataY[self.i]) )
            z=self.ioda.lev[I]
            time=self.ioda.tt[I]
            fcst=self.ioda.obs[I]-self.ioda.omf[I]
            ana=self.ioda.obs[I]-self.ioda.oma[I]
            obsi=self.ioda.obs[I]
            obsvarid=self.ioda.varid[I]
            obserrori=self.ioda.obserror[I]

            for uvarid in tqdm(np.unique(self.ioda.varid)):
                # Select variable
                II=np.where( obsvarid == uvarid )
                uz=z[II]
                ufcst=fcst[II]
                uana=ana[II]
                uobs=obsi[II]
                uobserror=obserrori[II]

                # Sort by depth
                III = np.argsort(uz)

                # Plot obs, background, analysis
                plot_prof(self.axis3, ufcst[III], uana[III], uobs[III], uz[III], var_name)

                # Plot omf, omaq
                plot_prof(self.axis4, uobs[III]-ufcst[III], uobs[III]-uana[III], ufcst[III]-uana[III], uz[III], var_name)

                # Add obs info to the figure
                self.axis5 = self.figure2.add_axes([0.75,0.75,0.2,0.2],frameon=False)
                self.axis5.axis('off')
                strtxt = '{0:10} {1}'.format('Instrument: ', inst_name) + '\n' + \
                         '{0:5} {1:3.2f}'.format('Lon:', self.dataX[self.i]) + '\n' + \
                         '{0:5} {1:3.2f}'.format('Lat:', self.dataY[self.i]) + '\n'
                self.axis5.text(0.01,0.55,strtxt,fontsize=20, fontweight='bold')

                self.fignum +=1

            plt.show()

        if event.button == 2:    # Middle mouse click: Regression plot for all instruments
            # Isolate variable type
            for INSTID in tqdm(np.unique(self.ioda.instid)):
                # Identify instrument
                for instrument in dict_inst:
                    if INSTID==dict_inst[instrument].instid:
                        inst_name=dict_inst[instrument].name

                I=np.where( (self.ioda.instid==INSTID) & (self.ioda.obs>=-99.0) & (-self.ioda.omf+self.ioda.obs!=0.0))
                figure2 = plt.figure(num=self.fignum, figsize=(16, 12), facecolor='c')

                axis2 = figure2.add_subplot(211)
                titlestr=inst_name+' OMF'
                plot_reg(axis2, self.ioda.obs[I], -self.ioda.omf[I]+self.ioda.obs[I], 'g', titlestr)

                axis3 = figure2.add_subplot(212)
                titlestr=inst_name+' OMA'
                plot_reg(axis3, self.ioda.obs[I], -self.ioda.oma[I]+self.ioda.obs[I], 'r', titlestr)

                self.fignum+=1

            plt.show()

        if event.button == 3:    # Right mouse click: Horizontal scatter plot of omf's and oma's for all instruments
            # Isolate var type
            for INSTID in tqdm(np.unique(self.ioda.instid)):
                # Identify instrument
                for instrument in dict_inst:
                    if INSTID==dict_inst[instrument].instid:
                        inst_name=dict_inst[instrument].name

                figure2 = plt.figure(num=self.fignum, figsize=(16, 12), facecolor='c')
                axis2 = figure2.add_subplot(211)
                map=draw_map()
                I=np.where(self.ioda.instid==INSTID)
                STD=np.std(self.ioda.omf[I])
                for shift in [0, 360]:
                    x, y =map(self.dataX[I]+shift,self.dataY[I])
                    axis2.scatter(x, y, 5, c=self.ioda.omf[I], cmap=cm.bwr,vmin=-2*STD,vmax=2*STD,edgecolor=None,lw=0)
                titlestr=inst_name+' OMF'
                plt.title(titlestr,fontsize=24,fontweight='bold')

                axis3 = figure2.add_subplot(212)
                map=draw_map()
                for shift in [0, 360]:
                    x, y =map(self.dataX[I]+shift,self.dataY[I])
                    axis3.scatter(x, y, 5, c=self.ioda.oma[I], cmap=cm.bwr,vmin=-2*STD,vmax=2*STD,edgecolor=None,lw=0)
                titlestr=inst_name+' OMA'
                plt.title(titlestr,fontsize=24,fontweight='bold')
                self.fignum+=1

            ax4 = figure2.add_axes([0.15, 0.25, 0.025, 0.5])
            norm = mpl.colors.Normalize(vmin=-2*STD,vmax=2*STD)
            mpl.colorbar.ColorbarBase(ax4, cmap=cm.bwr,norm=norm,orientation='vertical',extend='both')
            plt.show()

if __name__ == '__main__':
    description = """Observation space interactive map:
                     oceanview.py -i prof.out_*.nc adt.out_*.nc sst.out_*.nc"""

    parser = ArgumentParser(
        description=description,
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i',
        '--input',
        help="ioda files from the output of soca DA",
        type=str, nargs='+', required=True)

    print("""
             ============================================
             === Mouse left click: Profiles
             === Mouse middle click: regression
             === Mouse right click: Horizontal omf's/oma's
             === Usage: oceanview.py -i prof.out_*.nc adt.out_*.nc sst.out_*.nc
             ============================================
          """)
    args = parser.parse_args()
    listoffiles = args.input
    dict_inst, dict_varspecs, COLORS = dictionaries()
    example=observation_space(listoffiles)
    plt.show()
