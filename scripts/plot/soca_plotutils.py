from netCDF4 import Dataset, num2date, date2num
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as stats
plt.rcParams['axes.facecolor'] = 'xkcd:grey'
plt.rcParams['contour.negative_linestyle'] = 'solid'

def plothor(x,y,z,map,varname='',label='[m]',clim=[0,1],obs=[]):
    a=clim[0] #np.min(z[:])
    b=clim[1] #np.max(z[:])

    #print 'clim=',clim
    #print np.max(np.abs(z))

    clevs = np.linspace(a, b, 41)

    map.drawcoastlines()
    map.fillcontinents(color='coral')
    map.drawparallels(np.arange(-90.,120.,15.))
    map.drawmeridians(np.arange(0.,420.,30.))

    cmap=cm.hsv
    #cmap=cm.nipy_spectral    
    #cmap=cm.coolwarm
    #cmap=cm.seismic
    #cmap=cm.bwr
    #cmap=cm.Accent        
    cs = map.contourf(x, y, z, clevs, cmap = cmap, extend='both')
    #cs = map.pcolormesh(x, y, z, cmap = cmap)
    plt.title(varname)
    cbar=map.colorbar(shrink=0.2,format="%.1f")
    cbar.set_label(label, rotation=270)

    clevs = np.linspace(a, b, 21)
    map.contour(x, y, z, clevs, colors = 'k')
    
    xo, yo = map(obs.lon,obs.lat)
    map.plot(xo,yo,'.k',alpha=0.8,markersize=1)
    #map.scatter(xo,yo,c=obs.omb,alpha=0.9,edgecolors='none',cmap=cmap)
    #clevs = np.linspace(a, b, 11)
    #map.contour(x, y, z, clevs, colors='k')
    
class Grid:
    def __init__(self,fname='../../test/Data/360x210x63/ocean_geometry.nc'):
        ncfile = Dataset(fname,'r')        
        try:
            self.lat=np.squeeze(ncfile.variables['geolat'][:])
            self.lon=np.squeeze(ncfile.variables['geolon'][:])
            self.mask=np.squeeze(ncfile.variables['wet'][:])
        #print np.shape(self.wet)
        except:
            self.lat=np.squeeze(ncfile.variables['lat'][:])
            self.lon=np.squeeze(ncfile.variables['lon'][:])
            self.mask=np.squeeze(ncfile.variables['mask2d'][:])
        ncfile.close()        

class OceanState:
    def __init__(self, filename, maptype='N'):
        #print filename
        ncfile = Dataset(filename,'r')
        try:
            self.temp=np.squeeze(ncfile.variables['temp'][:])
            self.salt=np.squeeze(ncfile.variables['salt'][:])
            self.h=np.squeeze(ncfile.variables['h'][:])
            self.ssh=np.squeeze(ncfile.variables['ssh'][:])
            self.cicen=np.squeeze(ncfile.variables['cicen'][:])
            self.hicen=np.squeeze(ncfile.variables['hicen'][:])    
            ncfile.close()            
        except:
            self.temp=np.squeeze(ncfile.variables['Temp'][:])
            self.salt=np.squeeze(ncfile.variables['Salt'][:])
            self.ssh=np.squeeze(ncfile.variables['ave_ssh'][:])            
            self.h=np.squeeze(ncfile.variables['h'][:])
            
        self.grid=Grid()
        self.hbottomz=np.cumsum(self.h,axis=0)
        self.hmidz=self.hbottomz-0.5*self.h
        if (maptype=='N'):
            self.map = Basemap(projection='npstere',lon_0=0,boundinglat=50, resolution='l')
        elif (maptype=='S'):
            self.map = Basemap(projection='spstere',lon_0=0,boundinglat=50, resolution='l')
        else:            
            self.map = Basemap(projection='mill',lon_0=-100)
            #h = 30000.
            #self.map = Basemap(projection='nsper',lon_0=-35.0,lat_0=0,satellite_height=h*1000.,resolution='l')
            
        self.x, self.y = self.map(self.grid.lon,self.grid.lat)

        #self.obs=OceanObs(filename='/home/gvernier/Sandboxes/soca/soca-bundle-mom6/build-release/soca/test/fcst_da/Data/Jason-2-2018-04-15-soca-out_0000.nc')
        #self.obs.plot()

        self.obs=OceanObs(filename='/home/gvernier/Sandboxes/soca/soca-bundle-mom6/build-release/soca/test-ams/Data/Jason-2-2018-04-15-soca-out_0000.nc')
        #self.obs.plot()
        
    def plot_vert_section(self, other, fignum=1):
        j=100                         
        x=np.transpose(np.reshape(np.repeat(np.squeeze(self.grid.lon[j,:]),63),(360,63)))
        z=np.squeeze(-other.hmidz[:,j,:])        

        plt.figure(num=fignum)

        #plt.subplot(211)
        #print 'min temp:',np.min(self.temp[:,j,:]-other.temp[:,j,:])
        #print 'max temp:',np.max(self.temp[:,j,:]-other.temp[:,j,:])
                         
        incr = np.squeeze(self.temp[:,j,:])#-other.temp[:,j,:])

        depth=np.sum(other.h[:,j,:],axis=0)
        mask=np.ones(np.shape(depth))
        mask[np.where(depth<100.0)]=np.nan

        #print 'ocean levels: ',np.shape(incr)[0]
        for k in range(np.shape(incr)[0]):
            incr[k,:]=mask*incr[k,:]/self.grid.mask[j,:]

        #vmin=np.min(incr)
        #vmax=np.max(incr)
        vmin=0.0 #-3 #np.min(incr)
        vmax=31.0 #3 #abs(np.min(incr)) #np.max(incr)                
        clevs = np.linspace(vmin, vmax, 41)
        plt.contourf(x,z,incr, clevs, extend='both',cmap=cm.seismic)
        cbar=plt.colorbar(shrink=0.6,format="%.1f")
        
        clevs = np.linspace(vmin, vmax, 11)
        plt.contour(x,z,incr, clevs,colors='k')
        #plt.pcolor(x,z,incr,vmin=vmin,vmax=vmax,cmap=cm.seismic)
        #plt.pcolor(x-360,z,incr,vmin=vmin,vmax=vmax,cmap=cm.bwr)        
        plt.ylim((-6000, 0))
        #plt.xlim((-215, -195))

        #cbar.set_label('[K]', rotation=270)
        return
        plt.subplot(212)
        incr = self.salt[:,j,:]-other.salt[:,j,:]

        for k in range(np.shape(incr)[0]):
            incr[k,:]=mask*incr[k,:]/self.grid.mask[j,:]

        nk=np.shape(incr)[0]
        for iter in range(2):
            for k in np.arange(1,nk-1):
                incr[k,:]=(incr[k-1,:]+incr[k,:]+incr[k+1,:])/3.0

            
        #print "min salt incr:",np.min(incr)
        vmin=-0.2 
        vmax=0.2 #abs(np.min(incr)) #np.max(incr)        
        clevs = np.linspace(vmin, vmax, 41)
        #plt.pcolor(x,z,incr,vmin=vmin,vmax=vmax,cmap=cm.bwr)        
        plt.contourf(x,z,incr, clevs, extend='both',cmap=cm.seismic)
        
        #plt.pcolor(x,z,self.salt[:,j,:]-other.salt[:,j,:],vmin=-.2,vmax=.2,cmap=cm.bwr)
        plt.ylim((-6000, 0))
        #plt.xlim((-215, -195))
        cbar=plt.colorbar(shrink=0.6,format="%.1f")
        #cbar.set_label('[psu]', rotation=270)
        clevs = np.linspace(vmin, vmax, 11)
        plt.contour(x,z,incr, clevs,colors='k')

        
    def plot_horiz_section(self, other=[], vars=['temp'], levels=[0], fignum=1):
        plt.figure(num=fignum)
        #map = Basemap(projection='mill',lon_0=-100)
        #x, y = map(self.grid.lon,self.grid.lat)
        plothor(self.x,self.y,self.ssh,self.map,'',clim=[-1.6,1.4],label='',obs=self.obs)
        return
        for var in vars:
            if var=='temp':
                incr=self.temp[levels[0],:,:]-other.temp[levels[0],:,:]
                cmin=-5.5
                cmax=5.5
                titlestr='Temperature increment'
            if var=='salt':
                incr=self.salt[levels[0],:,:]-other.salt[levels[0],:,:]
                cmin=-.5
                cmax=.5
                titlestr='Salinity increment'
            if var=='ssh':
                incr=self.ssh[:,:]-other.ssh[:,:]
                cmin=-.3
                cmax=.3
                titlestr='SSH increment'
            if var=='cicen':
                self.maptype = 'N'
                incr=np.squeeze(np.sum(self.cicen[1:,:,:]-other.cicen[1:,:,:],0))
                cmin=-.5
                cmax=.5
                titlestr='cice increment'   
            if var=='hicen':
                self.maptype = 'N'
                incr=np.squeeze(np.sum(self.hicen[1:,:,:]-other.hicen[1:,:,:],0))
                cmin=-1.
                cmax=1.
                titlestr='hice increment'   
                
            plothor(self.x,self.y,incr,self.map,titlestr,clim=[cmin,cmax],label='',obs=self.obs)
            #plt.show()
        #plt.xlim((-215, -195))

    def plot_integrated_horiz(self, other, vars=['temp'], fignum=1):
        plt.figure(num=fignum)
        #map = Basemap(projection='mill',lon_0=-100)
        #x, y = map(self.grid.lon,self.grid.lat)

        for var in vars:
            if var=='temp':
                incr=np.sum(self.h[0:,:,:]*self.temp[0:,:,:],axis=0)/ \
                     np.sum(self.h[0:,:,:],axis=0) - \
                     np.sum(other.h[0:,:,:]*other.temp[0:,:,:],axis=0)/ \
                     np.sum(other.h[0:,:,:],axis=0)                     
                cmin=-0.02 #1300#np.min(incr)
                cmax=0.02 #1300#np.max(incr)
                titlestr='Temperature'
            if var=='salt':
                incr=np.sum(self.h[0:,:,:]*self.salt[0:,:,:],axis=0)/ \
                     np.sum(self.h[0:,:,:],axis=0) - \
                     np.sum(other.h[0:,:,:]*other.salt[0:,:,:],axis=0)/ \
                     np.sum(self.h[0:,:,:],axis=0)                     
                cmin=-0.01#np.min(incr)
                cmax=0.01 ##np.max(incr)
                titlestr='Salinity'                
                
            plothor(self.x,self.y,incr,self.map,titlestr,clim=[cmin,cmax],label='')

        
class OceanObs:
    def __init__(self, filename):

        ncfile = Dataset(filename,'r')
        self.lat=np.squeeze(ncfile.variables['latitude@GroupUndefined'][:])
        self.lon=np.squeeze(ncfile.variables['longitude@GroupUndefined'][:])    
        self.omb=-np.squeeze(ncfile.variables['obs_absolute_dynamic_topography@ombg'][:])
        self.oma=-np.squeeze(ncfile.variables['obs_absolute_dynamic_topography@oman'][:])
        self.obs=np.squeeze(ncfile.variables['obs_absolute_dynamic_topography@ObsValue'][:])
        ncfile.close()
        self.lon[self.lon>80.0]=self.lon[self.lon>80.0]-360

    def plot(self):
        print np.sum(np.abs(self.omb))/len(self.omb), stats.tstd(self.omb)
        print np.sum(np.abs(self.oma))/len(self.oma), stats.tstd(self.oma)


        map = Basemap(projection='mill',llcrnrlat=-72,urcrnrlat=72,
            llcrnrlon=-345,urcrnrlon=15,resolution='c')
#        map = Basemap(projection='mill',
#                      llcrnrlat = -70,
#                      urcrnrlat = 70)
        cmap=cm.nipy_spectral
        
        plt.subplot(121)        
        map.drawcoastlines()
        map.fillcontinents(color='coral')
        map.drawparallels(np.arange(-90.,120.,15.))
        map.drawmeridians(np.arange(0.,420.,30.))

        map.scatter(self.lon,self.lat,s=1,
                    c=self.omb,vmin=-0.3,vmax=0.3,cmap=cm.jet,latlon=True)

        plt.subplot(122)
        map.drawcoastlines()
        map.fillcontinents(color='coral')
        map.drawparallels(np.arange(-90.,120.,15.))
        map.drawmeridians(np.arange(0.,420.,30.))

        map.scatter(self.lon,self.lat,s=1,
                    c=self.oma,vmin=-0.3,vmax=0.3,cmap=cm.jet,latlon=True)
        #plt.colorbar()
        #plt.plot(self.oma,'.r',alpha=0.1)

        
        #mydensity = stats.gaussian_kde(data)        
        #X, Y = np.mgrid[-2:2:100j, -2:2:100j]
        #grid = np.vstack([X.ravel(), Y.ravel()])
        #Z = mydensity(grid)
        #Z = np.reshape(Z, X.shape)
        #plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth, extent=[-2,2,-2,2])
        #plt.plot(data[0], data[1], 'k.', markersize=2)
        #plt.xlim([-2,2])
        #plt.ylim([-2,2])
        
        plt.show()
