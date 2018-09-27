#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import glob
import matplotlib
from soca_plotutils import *

# Analysis increment
fname_ana='/home/gvernier/Sandboxes/soca/bmatrix2/soca-bundle/build-debug/soca/test/Data/3dvar.an.2018-04-15T00:00:00Z.nc'
fname_bkg='/home/gvernier/Sandboxes/soca/bmatrix2/soca-bundle/build-debug/soca/test/Data/example.fc.2018-04-15T00:00:00Z.P1D.nc'
oana=OceanState(fname_ana, maptype='R')
obkg=OceanState(fname_bkg)
oana.plot_horiz_section(obkg,vars=['temp','salt','ssh'],fignum=1)
plt.show()
oana.plot_vert_section(obkg,fignum=1)
plt.show()

#Reforecast increment
#fname_mom='MOM.res_Y2018_D109_S32400.nc'
#fname_ana='/home/gvernier/Sandboxes/soca/bmatrix2/soca-bundle/build/soca/test/scratch/RESTART/'+fname_mom
#fname_bkg='/home/gvernier/Sandboxes/soca/bmatrix2/soca-bundle/build/soca/test/scratch_noda/RESTART/'+fname_mom

#oana=OceanState(fname_ana)
#obkg=OceanState(fname_bkg)
#oana.plot_vert_section(obkg, fignum=2)

#plt.show()

#Obsspace
obs=OceanObs()
obs.plot_regress()
#obs.plot()
plt.show()
