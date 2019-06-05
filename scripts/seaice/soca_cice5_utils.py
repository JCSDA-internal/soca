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
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import scipy.io.netcdf as netcdf
import scipy
from matplotlib.colors import LinearSegmentedColormap
import soca_cice5_thermo as thermo

class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""
    
    name = 'nlcmap'
    
    def __init__(self, cmap, levels):
        self.cmap = cmap
        # @MRR: Need to add N for backend
        self.N = cmap.N
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels / self.levels.max()
        self._y = np.linspace(0.0, 1.0, len(self.levels))
    
    #@MRR Need to add **kw for 'bytes'
    def __call__(self, xi, alpha=1.0, **kw):
        """docstring for fname"""
        # @MRR: Appears broken? 
        # It appears something's wrong with the
        # dimensionality of a calculation intermediate
        #yi = stineman_interp(xi, self._x, self._y)
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi, alpha)

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
    print fname,varname
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

def plot2d(x, y, var, SUM=True, vmin=-0.1, vmax=0.1, CMAP=cm.jet, map=None, pole='N'):
    """Basic polar stereographic plot
    
    Args:
        x     (2D array float): Longitude array
        y     (2D array float): Latitude array
        var (2/3D array float): variable to be ploted
        
    Optional Args:
        SUM   (bool): Aggregate var
        vmin (float): min value for levels
        vmax (float): max value for levels
    """
    if map==None:
        #map = Basemap(projection='hammer',lon_0=180)
        #map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
        #map = Basemap(projection='ortho',lat_0=-45,lon_0=-100,resolution='l')
        if pole=='N':
            map = Basemap(projection='npstere',lon_0=0,boundinglat=50, resolution='l')
        if pole=='S':
            map = Basemap(projection='spstere',lon_0=0,boundinglat=-55, resolution='l')
        if pole=='C':        
            #map = Basemap(projection='mill',lon_0=180)
            map = Basemap(projection='hammer',lon_0=200)            
        map.drawcoastlines()
        map.fillcontinents(color='coral')
        map.drawparallels(np.arange(-90.,120.,15.))
        map.drawmeridians(np.arange(0.,420.,30.))
        
    #vmin=max(np.min(var.flatten()),-1e17)
    #vmax=min(np.max(var.flatten()),1e17)

    print 'vmin=',vmin
    print 'vmax=',vmax
    print '===================='
    #vmin = 0.5*vmin
    #vmin = 0.0 #-vmax
    #dl=(vmax-vmin)/25.0
    #levels=np.arange(vmin, vmax, dl)
    levels=np.linspace(vmin, vmax, 25, endpoint=True)

    #if SUM:
    #    VAR=np.sum(var, axis=0)
    #    VAR[np.abs(VAR)<0.0]=np.nan
    #    #VAR[VAR>1.0]=np.nan
    #else:
    VAR=var*(var/var)
    #map.contourf(x, y, VAR, levels, origin='lower',cmap=CMAP, extend='both',latlon=True)    
    try:
        map.contourf(x, y, VAR, levels, origin='lower',cmap=CMAP, extend='both',latlon=True)
        plt.colorbar(shrink=0.5)        
    except:
        pass
        
    #plt.imshow(VAR)
    return map

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
                print 'No mask'
            if ocnfname:
                print 'Reading ocean state ...'
                self.ocnfname=ocnfname
                self.Tocn = readvar(ocnfname, 'Temp'); self.Tocn=self.Tocn[0,:,:]
                self.Socn = readvar(ocnfname, 'Salt'); self.Socn=self.Socn[0,:,:]
                self.Tm = thermo.Tm(self.Socn)
        else:
            print 'No checkpoint file provided, initializing state to 0'
            self.ni=np.shape(self.x)[1]
            self.nj=np.shape(self.x)[0]
            self.ncat=5 # <--- hard coded, need to change that!
            print 'ni=',self.ni,' nj=',self.nj
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
        
    def __add__(self, other):
        spo=SeaIceState()
        spo.aicen[:] = self.aicen + other.aicen
        spo.vicen[:] = self.vicen + other.vicen
        spo.vsnon[:] = self.vsnon + other.vsnon
        spo.tsfcn[:] = self.tsfcn + other.tsfcn
        spo.hicen[:] = self.hicen + other.hicen
        spo.hsnon[:] = self.hsnon + other.hsnon
        spo.aice[:] = self.aice + other.aice
        spo.hice[:] = self.hice + other.hice
        spo.hsno[:] = self.hsno + other.hsno
        spo.tsfc[:] = self.tsfc + other.tsfc
        return spo

    def __sub__(self, other):
        spo=SeaIceState()
        spo.aicen[:] = self.aicen - other.aicen
        spo.vicen[:] = self.vicen - other.vicen
        spo.vsnon[:] = self.vsnon - other.vsnon
        spo.tsfcn[:] = self.tsfcn - other.tsfcn
        spo.hicen[:] = self.hicen - other.hicen
        spo.hsnon[:] = self.hsnon - other.hsnon
        spo.aice[:] = self.aice - other.aice
        spo.hice[:] = self.hice - other.hice
        spo.hsno[:] = self.hsno - other.hsno
        spo.tsfc[:] = self.tsfc - other.tsfc
        return spo
    
    def delta(self, ana, bkg):
        self.aicen[:] = ana.aicen - bkg.aicen
        self.vicen[:] = ana.vicen - bkg.vicen
        self.vsnon[:] = ana.vsnon - bkg.vsnon
        self.tsfcn[:] = ana.tsfcn - bkg.tsfcn
        self.hicen[:] = ana.hicen - bkg.hicen
        self.hsnon[:] = ana.hsnon - bkg.hsnon
        self.aice[:] = ana.aice - bkg.aice
        self.hice[:] = ana.hice - bkg.hice
        self.hsno[:] = ana.hsno - bkg.hsno
        self.tsfc[:] = ana.tsfc - bkg.tsfc

        print np.shape(self.aice)
        
        # case daice<0
        Im=np.where(self.aice<0.0)
        self.dT=np.zeros(np.shape(bkg.Tm))
        self.dT[Im]=self.aice[Im]*(-bkg.Tocn[Im]+bkg.Tm[Im])

        # case daice>=0
        alpha=1.0#
        Ip=np.where(self.aice>=0.0)
        self.dT[Ip]=alpha*self.aice[Ip]*(-bkg.Tocn[Ip]+bkg.Tm[Ip])
        print 'max dT=',np.max(self.dT), ' min dT=',np.min(self.dT)     

        #plt.imshow(self.dT,vmin=-0.15, vmax=0.15)
        #plt.show()
        #return self
        
    def info(self):
        print '================== ',self.descriptor,' ==================='
        for n in range(self.ncat):
            print 'cat=',str(n+1),' hi_max=',np.max(self.hicen[n,:,:].flatten()), \
                                  ' hi_min=',np.min(self.hicen[n,:,:].flatten())
            print '      hsno_max=',np.max(self.hsnon[n,:,:].flatten()), \
                                  ' hsno_min=',  np.min(self.hsnon[n,:,:].flatten())

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
        print 'filtering vsnon:',Imin
        self.vsnon[Imin]=0.0

        # We don't want (know how!) to create new ice ==> Mask out new ice
        if mask:
            for cat in range(self.ncat):
                self.aicen[cat,:,:]=self.aicen[cat,:,:]*self.iceumask
                self.vicen[cat,:,:]=self.vicen[cat,:,:]*self.iceumask
                self.vsnon[cat,:,:]=self.vsnon[cat,:,:]*self.iceumask

        # Make sure aggregate of aicen is in [0,1]        
        Imax=np.where(np.sum(self.aicen, axis=0)>1.)    
        self.aicen[:,Imax[0],Imax[1]]=self.aicen[:,Imax[0],Imax[1]]/np.sum(self.aicen[:,Imax[0],Imax[1]], axis=0)
        
        # By definition, surface temperature of snow or ice is in ]-alot, 0]
        self.tsfcn[self.tsfcn>0.0]=0.0

        # Update thicknesses
        I=np.where(self.aicen>1e-3)
        self.hicen[I] = self.vicen[I]/self.aicen[I]
        self.hsnon = np.zeros(np.shape(self.aicen))
        self.hsnon[I] = self.vsnon[I]/self.aicen[I]

        # Update aggregates
        self.aice=np.sum(self.aicen, axis=0)
        self.vice=np.sum(self.vicen, axis=0)
        self.hice=self.vice/self.aice
        self.hsno=self.vsno/self.aice
        self.tsfc=np.sum(self.aicen*self.tsfcn, axis=0)

    def update_ocn(self):
        # Assert we operate on increment
        self.dT=self.aice*(self.Tm-self.To)
                
    def hist(self, varname='aicen', type='incr', CMAP=cm.BrBG):

        for n in range(self.ncat):
            plt.figure(num=n+1)
            I=np.where( (self.y>55.0) & (self.aicen[n,:,:]>0.0) )
            #plt.plot(scipy.special.logit(self.aicen[n,I[0],I[1]].flatten()),\
            #         scipy.special.logit(self.hicen[n,I[0],I[1]].flatten()),'.g', alpha=0.1)
            nn, bins, patches = plt.hist(scipy.special.logit(self.aicen[n,I[0],I[1]].flatten()),\
                                     500, normed=1, facecolor='green', alpha=0.75)
        plt.show()

    def plot(self, varname='aicen', type='incr', CMAP=cm.jet, NLMAP=False, vmin=-0.1, vmax=0.1, figname='fig.png'):
        if varname=='aicen': var=self.aicen; agg=self.aice
        if varname=='hicen': var=self.hicen; agg=self.hice
        if varname=='hsnon': var=self.hsnon; agg=self.hsno
        if varname=='Tsfcn': var=self.tsfcn; agg=self.tsfc
        #if varname=='dT': var=self.dT; agg=self.tsfc
        
        plt.figure(num=1,figsize=(18,8))
        for n in range(self.ncat):
            print 'cat:',n+1,' ',np.shape(var),np.shape(self.x)
            plt.subplot(2,3,n+1)
            #vmin=np.min(var[n,:,:].flatten())
            #vmax=np.max(var[n,:,:].flatten())
            plot2d(self.x, self.y, var[n,:,:], vmin=vmin, vmax=vmax, CMAP=CMAP)
            #plt.colorbar(shrink=0.5)
            plt.title('cat-'+str(n+1))
        plt.subplot(2,3,n+2)
        if NLMAP:
            levels0=np.linspace(vmin, vmax, num=9, endpoint=True)
            print levels0
            amp=vmax-vmin
            levels=amp*np.tanh((levels0-0.5*amp)/amp)
            print levels
            cmap_nonlin = nlcmap(CMAP, levels)
            CMAP = cmap_nonlin
        print 'agg:',np.shape(agg),np.shape(self.x)
        #plot2d(self.x, self.y, agg, vmin=vmin, vmax=vmax, CMAP=CMAP)
        plot2d(self.x, self.y, self.dT, vmin=-0.15, vmax=0.15, CMAP=CMAP)
        #plt.colorbar(shrink=0.5)
        #plt.colorbar(spacing='proportional',
        #             orientation='vertical',
        #             shrink=0.7)#, format="%.0f",fontsize=28)
        #cbar.set_label(r"ET [mm/month]", size=10)



        #plt.title('Aggregate')
        plt.title('dT')
        plt.suptitle(varname)
        plt.savefig(figname)
        #plt.show()
        plt.clf()




