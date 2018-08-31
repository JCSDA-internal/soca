from netCDF4 import Dataset, num2date, date2num
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class Grid:
    def __init__(self):    
        fname='/home/gvernier/Sandboxes/soca/soca-bundle/soca/test/Data/360x210x63/ocean_geometry.nc'
        ncfile = Dataset(fname,'r')
        self.lat=np.squeeze(ncfile.variables['geolat'][:])
        self.lon=np.squeeze(ncfile.variables['geolon'][:])         
        ncfile.close()        

class OceanState:
    def __init__(self, filename):
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
            self.h=np.squeeze(ncfile.variables['h'][:])
            
        self.grid=Grid()
        self.hbottomz=np.cumsum(self.h,axis=0)
        self.hmidz=self.hbottomz-0.5*self.h        
        self.map = Basemap(projection='mill',lon_0=-100)
        self.x, self.y = self.map(self.grid.lon,self.grid.lat)

    def plot_vert_section(self, other, fignum=1):
        x=np.transpose(np.reshape(np.repeat(np.squeeze(self.grid.lon[100,:]),63),(360,63)))
        #z=np.squeeze(-self.hmidz[:,100,:])
        z=np.squeeze(-other.hmidz[:,100,:])        

        plt.figure(num=fignum)
        j=90
        plt.subplot(211)
        plt.pcolor(x,z,self.temp[:,j,:]-other.temp[:,j,:],vmin=-5.0,vmax=5.0,cmap=cm.bwr)
        plt.ylim((-1000, 0))
        #plt.xlim((-215, -195))

        plt.subplot(212)
        plt.pcolor(x,z,self.salt[:,j,:]-other.salt[:,j,:],vmin=-.2,vmax=.2,cmap=cm.bwr)
        plt.ylim((-1000, 0))
        #plt.xlim((-215, -195))

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
        
