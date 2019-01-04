#! /usr/bin/env python

from netCDF4 import Dataset, num2date, date2num
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import glob
import matplotlib

map = Basemap(projection='mill',lon_0=-100)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
matplotlib.rc('font', **font)

fname='../../../build/soca/test/geom_output.nc'
ncfile = Dataset(fname,'r')
x=np.squeeze(ncfile.variables['lon'][:])
x[x<0]=x[x<0]+360
y=np.squeeze(ncfile.variables['lat'][:])
mask=np.squeeze(ncfile.variables['mask2d'][:])
ncfile.close()

x=np.reshape(x,75600)
y=np.reshape(y,75600)
mask=np.reshape(mask,75600)
x,y=map(x,y)
I=np.where(mask==1)
triang = tri.Triangulation(x[I], y[I])
#triang = tri.Triangulation(x, y, mask)

#mask = np.reshape(mask,75600)
#print np.shape(x), np.shape(y), np.shape(mask)
#triang.set_mask(mask)

#plt.triplot(triang, 'bo-', lw=1)

#plt.triplot(triang, '-k', lw=1)

#map.drawcoastlines()
#map.fillcontinents(color='coral')
#map.drawparallels(np.arange(-90.,120.,15.))
#map.drawmeridians(np.arange(0.,420.,30.))


#plt.show()

flist=glob.glob('../../../build-release/soca/test/geom_output_*.nc')
flist.sort()

obsflist=glob.glob('../../../build/soca/test/adt-test.nc_*.nc')
#obsflist=glob.glob('../../../build/soca/test/sic-test.nc_*.nc')
obsflist.sort()


map.drawcoastlines()
#map.fillcontinents(color='coral')
map.drawparallels(np.arange(-90.,120.,15.))
map.drawmeridians(np.arange(0.,420.,30.))

color=['r','k','b','y','g','m','c','y','r','k','b','y','g','m','c','y']
marker=['+','.','o','*','o','.','o','.','o','.','o','.''o','.','o','.']
cnt=0
for fname in flist:
    if cnt>15:
        cnt=0
    # 2d pe domain
    print fname
    ncfile = Dataset(fname,'r')
    x=np.squeeze(ncfile.variables['lon'][:])
    x[x<0]=x[x<0]+360
    y=np.squeeze(ncfile.variables['lat'][:])
    mask=np.squeeze(ncfile.variables['mask'][:])
    shoremask=np.squeeze(ncfile.variables['shoremask'][:])    
    ncfile.close()
    I=np.where(mask==1)
    #Is=np.where(abs(shoremask-mask)!=0)
    Is=np.where(mask==0)    
    #Is=np.where(shoremask==1)    
    #plt.subplot(121)
    map.plot(x[I],y[I],color='k',marker='.',linestyle='None',alpha=0.8,latlon=True)    
    #map.plot(x[I],y[I],color=color[cnt],marker='.',linestyle='None',alpha=0.8,latlon=True)
    #plt.subplot(122)    
    #map.plot(x[Is],y[Is],color=color[cnt],marker='*',linestyle='None',alpha=0.1,latlon=True)
    #plt.plot(x[Is],y[Is],color=color[cnt],marker='.',linestyle='None',alpha=0.1)    
    try:
        # obs in pe domain
        print obsflist[cnt]
        ncfile = Dataset(obsflist[cnt],'r')
        xo=np.squeeze(ncfile.variables['lon'][:])
        #xo[xo<0]=xo[xo<0]+360    
        yo=np.squeeze(ncfile.variables['lat'][:])
        ncfile.close()    

        #plt.subplot(122)
        map.plot(xo,yo,color[cnt],marker=marker[cnt],linestyle='None',alpha=0.5,latlon=True)
    except:
        print 'oops ...'
    cnt+=1

plt.show()

