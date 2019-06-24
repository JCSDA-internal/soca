# -*- coding: utf-8 -*-
"""CICE5 utilies and class definition

This module contains classes, reading and writing function relevant to CICE5 

TODO:
     * Consolidate write and write2d
     * Document!
     * Add option for S polar stereo
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.io.netcdf as netcdf
import scipy

def readvar(fname, varname):
    """Function reading netcdf file
    
    Args:
        fname   (str): Full name (including path) of the netcdf file
        varname (str): Variable to be read from fname

    Return:
        If reading CICE5 checkpoint file:
            VAR (float), ni(int), nj(int), ncat(int)
        else:
            VAR (float)
    """
    ncfile = Dataset(fname,'r')
    VAR=np.squeeze(ncfile.variables[varname][:])
    try:
        ni=ncfile.dimensions['ni'].size
        nj=ncfile.dimensions['nj'].size
        ncat=ncfile.dimensions['ncat'].size   
    except:
        pass
    ncfile.close()
    try:
        return VAR, ni, nj, ncat
    except:
        return VAR

class SeaIceState:
    def __init__(self, gridname='geom_output.nc', fname=None, ocnfname=None, descriptor='Sea-ice state'):
        self.fname = fname
        self.x = readvar(gridname, 'lon')
        self.y = readvar(gridname, 'lat')
        self.descriptor=descriptor

        # Read or initialize sea-ice state
        if fname:
            self.aicen, self.ni, self.nj, self.ncat = readvar(fname, 'aicen')
            self.vicen, self.ni, self.nj, self.ncat = readvar(fname, 'vicen')
            self.vsnon, self.ni, self.nj, self.ncat = readvar(fname, 'vsnon')
            self.tsfcn, self.ni, self.nj, self.ncat = readvar(fname, 'Tsfcn')
            try:
                self.iceumask, self.ni, self.nj, self.ncat = readvar(fname, 'iceumask')
            except:
                print('No mask')
            if ocnfname:
                print('Reading ocean state ...')
                self.ocnfname=ocnfname
                self.Tocn = readvar(ocnfname, 'Temp'); self.Tocn=self.Tocn[0,:,:]
                self.Socn = readvar(ocnfname, 'Salt'); self.Socn=self.Socn[0,:,:]
                self.Tm = -0.054*self.Socn
        else:
            print('No checkpoint file provided, initializing state to 0')
            self.ni=np.shape(self.x)[1]
            self.nj=np.shape(self.x)[0]
            self.ncat=5 # <--- hard coded, need to change that!
            print('ni=',self.ni,' nj=',self.nj)
            shape=(self.ncat,self.nj,self.ni)
            self.aicen=np.zeros(shape)
            self.vicen=np.zeros(shape)
            self.vsnon=np.zeros(shape)
            self.tsfcn=np.zeros(shape)
            self.iceumask=np.zeros(np.shape(self.x.T))

        # Compute diagnositcs from sea-ice state        
        self.hicen = np.zeros(np.shape(self.aicen))
        I=np.where(self.aicen>0.0)
        self.hicen[I] = self.vicen[I]/self.aicen[I]
        self.hsnon = np.zeros(np.shape(self.aicen))
        self.hsnon[I] = self.vsnon[I]/self.aicen[I]

        # Compute aggregates
        self.aice=np.sum(self.aicen, axis=0)
        self.vice=np.sum(self.vicen, axis=0)
        self.vsno=np.sum(self.vsnon, axis=0)
        self.hice=self.vice/self.aice
        self.hsno=self.vsno/self.aice
        self.tsfc=np.sum(self.aicen*self.tsfcn, axis=0)
        self.iceumask=np.zeros(np.shape(self.aice))
        self.iceumask[self.aice>0.0]=1.0

    def info(self):
        print('================== ',self.descriptor,' ===================')
        print('aice_max=',np.max(np.max(np.max(np.sum(self.aicen[:,:,:],axis=0)))))
        print('tsfcn_sum=',np.sum(self.tsfcn.flatten()))
        for n in range(self.ncat):
            print('cat=',str(n+1),' hi_max=',np.max(self.hicen[n,:,:].flatten()), \
                                  ' hi_min=',np.min(self.hicen[n,:,:].flatten()))
            print('      hsno_max=',np.max(self.hsnon[n,:,:].flatten()), \
                                  ' hsno_min=',  np.min(self.hsnon[n,:,:].flatten()))

    def write(self, VAR, fname, varname, mode='a', size='3D'):        
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

    def write2d(self, VAR, fname, varname, mode='a'):
        ncfile = netcdf.netcdf_file(fname,mode)
        VARtmp=ncfile.variables[varname][:,:]
        VARtmp[:]=VAR[:]
        ncfile.close()
        
    def update(self):
        self.aice=np.sum(self.aicen, axis=0)
        self.vice=np.sum(self.vicen, axis=0)
        self.hicen=self.vicen/self.aicen
        self.hsnon=self.vicen/self.aicen
        self.hice=self.vice/self.aice
        self.hsno=self.vsno/self.aice
        self.tsfc=np.sum(self.aicen*self.tsfcn, axis=0)

    def fix_bounds(self,mask=False):
        # aicen is in [0,1]
        vmin=0.0
        vmax=1.0
        Imin=np.where(self.aicen<vmin)
        Imax=np.where(self.aicen>vmax)        
        self.aicen[Imin]=vmin
        self.aicen[Imax]=vmax

        # vicen>=0
        Imin=np.where(self.vicen<0.0)
        self.vicen[Imin]=0.0

        # vsnon>=0
        Imin=np.where(self.vsnon<0.0)
        self.vsnon[Imin]=0.0

        # Make sure aggregate of aicen is in [0,1]        
        Imax=np.where(np.sum(self.aicen, axis=0)>1.0)    
        self.aicen[:,Imax[0],Imax[1]]=self.aicen[:,Imax[0],Imax[1]]/(np.sum(self.aicen[:,Imax[0],Imax[1]], axis=0)+1e-8)
        
        # By definition, surface temperature of snow or ice is in ]-alot, 0]
        self.tsfcn[self.tsfcn>0.0]=0.0

        # Update thicknesses
        #I=np.where(self.aicen>1e-3)
        #self.hicen[I] = self.vicen[I]/self.aicen[I]
        #self.hsnon = np.zeros(np.shape(self.aicen))
        #self.hsnon[I] = self.vsnon[I]/self.aicen[I]

        # Update aggregates
        self.aice=np.sum(self.aicen, axis=0)
        self.vice=np.sum(self.vicen, axis=0)
        self.hice=self.vice/self.aice
        self.hsno=self.vsno/self.aice
        self.tsfc=np.sum(self.aicen*self.tsfcn, axis=0)

    def update_ocn(self):
        # Assert we operate on increment
        self.dT=self.aice*(self.Tm-self.To)
