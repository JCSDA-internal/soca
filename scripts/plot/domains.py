#! /usr/bin/env python

from netCDF4 import Dataset, num2date, date2num
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import glob
import matplotlib


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
matplotlib.rc('font', **font)

flist=glob.glob('../../../build/soca/test/geom_output_*.nc')
flist.sort()

obsflist=glob.glob('../../../build/soca/test/adt-test.nc_*.nc')
#obsflist=glob.glob('../../../build/soca/test/sic-test.nc_*.nc')
obsflist.sort()

color=['r','k','b','y','g','m','c','y','r','k','b','y','g','m','c','y']
marker=['+','.','o','*','o','.','o','.','o','.','o','.''o','.','o','.']
cnt=0
for fname in flist:
    if cnt>15:
        cnt=0
    # 2d pe domain
    ncfile = Dataset(fname,'r')
    x=np.squeeze(ncfile.variables['lon'][:])
    x[x<0]=x[x<0]+360
    y=np.squeeze(ncfile.variables['lat'][:])
    mask=np.squeeze(ncfile.variables['obsmask'][:])
    ncfile.close()
    I=np.where(mask==1)
    #plt.subplot(121)    
    plt.plot(x[I],y[I],color=color[cnt],marker='.',linestyle='None',alpha=0.1)

    try:
        # obs in pe domain
        print obsflist[cnt]
        ncfile = Dataset(obsflist[cnt],'r')
        xo=np.squeeze(ncfile.variables['lon'][:])
        #xo[xo<0]=xo[xo<0]+360    
        yo=np.squeeze(ncfile.variables['lat'][:])
        ncfile.close()    

        #plt.subplot(122)
        plt.plot(xo,yo,color[cnt],marker=marker[cnt],linestyle='None',alpha=0.5)
    except:
        print 'oops ...'
    cnt+=1
    
plt.show()
