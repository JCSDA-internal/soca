from netCDF4 import Dataset, num2date, date2num
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plothor(x,y,z,map,varname='',label='[m]',clim=[0,1]):
    a=clim[0] #np.min(z[:])
    b=clim[1] #np.max(z[:])

    print 'clim=',clim
    print np.max(np.abs(z))

    clevs = np.linspace(a, b, 41)

    map.drawcoastlines()
    map.fillcontinents(color='coral')
    map.drawparallels(np.arange(-90.,120.,15.))
    map.drawmeridians(np.arange(0.,420.,30.))
    #cmap=cm.nipy_spectral
    cmap=cm.coolwarm
    #cmap=cm.seismic
    #cmap=cm.bwr
    #cmap=cm.Accent        
    cs = map.contourf(x, y, z, clevs, cmap = cmap, extend='both')
    #cs = map.pcolormesh(x, y, z, cmap = cmap)
    plt.title(varname)
    cbar=plt.colorbar(shrink=0.5,format="%.1f")
    cbar.set_label(label, rotation=270)
    #clevs = np.linspace(a, b, 11)
    #map.contour(x, y, z, clevs, colors='k')
    
class Grid:
    def __init__(self):    
        fname='../../test/Data/360x210x63/ocean_geometry.nc'
        ncfile = Dataset(fname,'r')
        self.lat=np.squeeze(ncfile.variables['geolat'][:])
        self.lon=np.squeeze(ncfile.variables['geolon'][:])         
        ncfile.close()        

class OceanState:
    def __init__(self, filename, maptype='N'):
        print filename
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
        self.x, self.y = self.map(self.grid.lon,self.grid.lat)

    def plot_vert_section(self, other, fignum=1):
        x=np.transpose(np.reshape(np.repeat(np.squeeze(self.grid.lon[100,:]),63),(360,63)))
        #z=np.squeeze(-self.hmidz[:,100,:])
        z=np.squeeze(-other.hmidz[:,100,:])        

        plt.figure(num=fignum)
        j=110
        plt.subplot(211)
        print 'min temp:',np.min(self.temp[:,j,:]-other.temp[:,j,:])
        print 'max temp:',np.max(self.temp[:,j,:]-other.temp[:,j,:])
        incr = self.temp[:,j,:]-other.temp[:,j,:]
        #vmin=np.min(incr)
        #vmax=np.max(incr)
        vmin=-0.1 #np.min(incr)
        vmax=0.1 #abs(np.min(incr)) #np.max(incr)                
        clevs = np.linspace(vmin, vmax, 41)
        plt.contourf(x,z,incr, clevs, extend='both',cmap=cm.spectral)
        #plt.pcolor(x,z,self.temp[:,j,:]-other.temp[:,j,:],vmin=-.05,vmax=.05,cmap=cm.bwr)
        plt.ylim((-2000, 0))
        #plt.xlim((-215, -195))
        cbar=plt.colorbar(shrink=0.5,format="%.1f")
        cbar.set_label('[K]', rotation=270)

        plt.subplot(212)
        incr = self.salt[:,j,:]-other.salt[:,j,:]        
        vmin=-0.1 #np.min(incr)
        vmax=0.1 #abs(np.min(incr)) #np.max(incr)        
        clevs = np.linspace(vmin, vmax, 41)        
        plt.contourf(x,z,incr, clevs, extend='both',cmap=cm.spectral)
        
        #plt.pcolor(x,z,self.salt[:,j,:]-other.salt[:,j,:],vmin=-.2,vmax=.2,cmap=cm.bwr)
        plt.ylim((-2000, 0))
        #plt.xlim((-215, -195))
        cbar=plt.colorbar(shrink=0.5,format="%.1f")
        cbar.set_label('[psu]', rotation=270)
        
    def plot_horiz_section(self, other, vars=['temp'], levels=[0], fignum=1):
        plt.figure(num=fignum)
        #map = Basemap(projection='mill',lon_0=-100)
        #x, y = map(self.grid.lon,self.grid.lat)

        for var in vars:
            if var=='temp':
                incr=self.temp[levels[0],:,:]-other.temp[levels[0],:,:]
                cmin=-0.1#2.5
                cmax=0.1#2.5
                titlestr='Temperature increment'
            if var=='salt':
                incr=self.salt[levels[0],:,:]-other.salt[levels[0],:,:]
                cmin=-.2
                cmax=.2
                titlestr='Salinity increment'
            if var=='ssh':
                incr=self.ssh[:,:]-other.ssh[:,:]
                cmin=-.2
                cmax=.2
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
                
            plothor(self.x,self.y,incr,self.map,titlestr,clim=[cmin,cmax],label='')
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
                cmin=-0.05 #1300#np.min(incr)
                cmax=0.05 #1300#np.max(incr)
                titlestr='Temperature'
            if var=='salt':
                incr=np.sum(self.h[0:,:,:]*self.salt[0:,:,:],axis=0)/ \
                     np.sum(self.h[0:,:,:],axis=0) - \
                     np.sum(other.h[0:,:,:]*other.salt[0:,:,:],axis=0)/ \
                     np.sum(self.h[0:,:,:],axis=0)                     
                cmin=-0.05#np.min(incr)
                cmax=0.05 ##np.max(incr)
                titlestr='Salinity'                
                
            plothor(self.x,self.y,incr,self.map,titlestr,clim=[cmin,cmax],label='')

        
class OceanObs:
    def __init__(self, basedir='/home/gvernier/Sandboxes/soca/bmatrix2/soca-bundle/build/soca/test/',basename='40'):

        # Observations
        fname=basedir+'fort.'+basename+'1'
        f = open(fname, 'r')
        cnt=0
        for line in f:
            cnt+=1
        f.close()
        self.nobs=cnt
        f = open(fname, 'r')
        self.obs=np.zeros(cnt)
        self.lon=np.zeros(cnt)
        self.lat=np.zeros(cnt)
        self.depth=np.zeros(cnt)
        cnt=0
        for line in f:
            line = line.strip()
            columns = line.split()
            self.lon[cnt] = float(columns[0])
            self.lat[cnt] = float(columns[1])
            self.obs[cnt] = float(columns[2]) #/100.0
            self.depth[cnt] = float(columns[3]) #/100.0    
            cnt+=1
        f.close()
        self.obs[self.obs<=-900.0]=np.nan
        # Simulated Observations
        fname=basedir+'fort.'+basename+'2'
        f = open(fname, 'r')
        cnt=0
        for line in f:
            cnt+=1
        f.close()

        f = open(fname, 'r')
        self.model=np.zeros(cnt)
        cnt=0
        for line in f:
            line = line.strip()
            columns = line.split()
            self.model[cnt] = float(columns[0])
            cnt+=1
        f.close()
        
    def plot(self):
        nobs=self.nobs
        nouter=len(self.model)/nobs
        print 'nouter',nouter
        print 'nobs=',self.nobs
        plt.plot(self.obs,-self.depth,'*b')
        plt.plot(self.model[0:nobs],-self.depth,'.-g')
        for iter in range(nouter-1):
            print 'bkg=',self.model[(iter+1)*nobs:(iter+2)*nobs]

            plt.plot(self.model[(iter+1)*nobs:(iter+2)*nobs],-self.depth,'.-r',lw=0.1)
        #plt.plot(bkg[:,1:],'-r')
        plt.grid(True)
        plt.show()

    def plot_regress(self):
        nobs=self.nobs
        nouter=len(self.model)/nobs
        print 'nouter',nouter
        print 'nobs=',self.nobs

        #plt.plot(self.model[0:nobs],self.obs,'.k')
        #plt.grid(True)
        #plt.show()

        colors=['g','r','m','y','b']
        for iter in range(nouter-1):
            plt.plot(self.model[(iter+1)*nobs:(iter+2)*nobs],self.obs,'.',alpha=0.1)#,colors=colors[iter])
        plt.grid(True)
        plt.show()
        
