#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute standard deviation of ensemble

This script reads an ensemble of states and computes the standard deviation of the
ensemble.

Example:
     $  ./_.................

Todo:

"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import scipy.io.netcdf as netcdf
import scipy
import sys
import glob
from datetime import datetime, timedelta
import matplotlib
import scipy.special as scsp

class State:
    def __init__(self, fname, varname):
        self.fname=fname
        self.varname=varname

    def read(self):
        print self.varname
        print self.fname        
        ncfile = Dataset(self.fname,'r')
        self.var=np.squeeze(ncfile.variables[self.varname][:])
        ncfile.close()
        self.ni=np.shape(self.var)[1]
        self.nj=np.shape(self.var)[0]

        
    def write(self, fname, mode='a', size='3D'):        
        ncfile = netcdf.netcdf_file(fname,mode)
        ni=self.ni
        nj=self.nj
        
        ncat=self.ncat
        if (mode!='a'):
            ncfile.createDimension('ni',ni)
            ncfile.createDimension('nj',nj)
            ncfile.createDimension('ncat',ncat)
        if (size=='3D'):
            VARtmp = ncfile.createVariable(varname,np.dtype('float32').char,('ncat','nj','ni'))
            VARtmp=ncfile.variables[varname][:,:,:]            
        if (size=='2D'):
            VARtmp = ncfile.createVariable(varname,np.dtype('float32').char,('nj','ni'))
            VARtmp=ncfile.variables[varname][:,:]

        VARtmp[:]=VAR[:]
        ncfile.close()
        

        
class EnsState:
    def __init__(self, flist,varname):
        self.flist = glob.glob(flist)
        self.varname = varname
        self.mean_initialized = False
        
    def ensemble_mean(self):
        cnt = 0.0
        for fname in self.flist:
            if cnt == 0.0:
                mean= State(fname, self.varname)
                mean.read()
            else:
                tmp = State(fname, self.varname)
                tmp.read()
                mean.var = mean.var + tmp.var
            cnt+=1
        self.mean = mean.var/cnt
        self.init_mean = True

    def ensemble_std(self):
        if (not(self.mean_initialized)):
            self.ensemble_mean()
        cnt = 0.0
        for fname in self.flist:
            if cnt == 0.0:
                std = State(fname, self.varname)
                std.read()
                std.var = (std.var - self.mean)**2
            else:
                std = State(fname, self.varname)
                std.read()
                std.var = std.var + (std.var - self.mean)**2

            cnt+=1
        self.std = std.var/cnt
        self.init_mean = True

        try:
            plt.pcolor(self.std[0,:,:],vmin=0,vmax=1.5)
        except:
            plt.pcolor(self.std[:,:],vmin=0,vmax=.05)            
        plt.colorbar()
        plt.show()
        
if __name__ == '__main__':

    print sys.argv[0],': Moments of ensemble'
    try:
        filelist=sys.argv[1]      # path to file with wild card. ex: ./ens/ens-0*.nc
    except:
        filelist = '/home/gvernier/Sandboxes/MOM6-examples/ice_ocean_SIS2/SIS2/scratch_ens/RESTART-??/MOM.res.nc'

    ssh_ens = EnsState(flist=filelist, varname='ave_ssh')
    ssh_ens.ensemble_std()        
    T_ens = EnsState(flist=filelist, varname='Temp')
    T_ens.ensemble_std()
    S_ens = EnsState(flist=filelist, varname='Salt')
    S_ens.ensemble_std()    
    
'''    cnt = 0.0
    for fname in flist:
        print fname
        if cnt == 0.0:
            ssh
        ncfile = Dataset(fname,'r')
        ssh_bkg=np.squeeze(ncfile.variables['ave_ssh'][:])
        t_bkg=np.squeeze(ncfile.variables['Temp'][:])
        s_bkg=np.squeeze(ncfile.variables['Salt'][:])
        ncfile.close()
    else:
        ncfile = Dataset(fname,'r')        
        ssh_bkg=ssh_bkg + np.squeeze(ncfile.variables['ave_ssh'][:])
        t_bkg=t_bkg + np.squeeze(ncfile.variables['Temp'][:])
        s_bkg=s_bkg + np.squeeze(ncfile.variables['Salt'][:])
        ncfile.close()        
    cnt+=1.0

ssh_bkg=ssh_bkg/cnt
t_bkg=t_bkg/cnt
s_bkg=s_bkg/cnt

cnt = 0.0
for fname in flist:
    print fname
    if cnt == 0.0:
        ncfile = Dataset(fname,'r')
        std_ssh=(np.squeeze(ncfile.variables['ave_ssh'][:])-ssh_bkg)**2
        #std_t=np.squeeze(ncfile.variables['Temp'][:])
        #std_s=np.squeeze(ncfile.variables['Salt'][:])
        ncfile.close()
    else:
        ncfile = Dataset(fname,'r')        
        std_ssh=std_ssh+(np.squeeze(ncfile.variables['ave_ssh'][:])-ssh_bkg)**2
        #std_t=t_bkg + np.squeeze(ncfile.variables['Temp'][:])
        #std_s=s_bkg + np.squeeze(ncfile.variables['Salt'][:])
        ncfile.close()        
    cnt+=1.0
    
std_ssh=std_ssh/(cnt-1.0)

plt.pcolor(std_ssh,vmin=0,vmax=0.05)
plt.colorbar()
plt.show()
'''
