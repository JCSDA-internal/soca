#!/usr/bin/env python

#
# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

"""CICE5 Inititalization

This script reads an analysis and its corresponding background 
and writes a properly balanced restart file. 

Example:
     $  ./soca_cice5_init.py gs04020.cice5_model.res.nc anal020.cice5_model.res.nc ogrid05.nc TEST.res.nc

TODO:
     * Currently not using the vicen increment
     * Modify mixed layer temperature accordingly to be consistent with the updated ice fraction
"""

from __future__ import print_function
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
import shutil
import soca_cice5_utils
import matplotlib.cm as cm
import warnings
import scipy.io.netcdf as netcdf
warnings.simplefilter("ignore", RuntimeWarning)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(
        description=('')
    )

    required = parser.add_argument_group(title='required arguments')
    required.add_argument(
        '-bkg', '--bkg_fname',
        help="CICE5 background restart file",
        type=str, nargs='+', required=True)
    required.add_argument(
        '-ocn', '--ocn_fname',
        help="MOM6 background or analysis restart file",
        type=str, nargs='+', required=True)    
    required.add_argument(
        '-ana', '--ana_fname',
        help="Sea-ice analysis file for aicen and hicen in the same grid as the background restart",
        type=str, nargs='+', required=True)
    required.add_argument(
        '-grid', '--grid_fname',
        help="Common grid for ocean and sea-ice",
        type=str, nargs='+', required=True)
    required.add_argument(
        '-o', '--output',
        help="Filename for the CICE5 restart containing the analysis",
        type=str, required=True)

    args = parser.parse_args()
    bkg_fname    = args.bkg_fname[0]
    ocn_fname    = args.ocn_fname[0] 
    ana_fname    = args.ana_fname[0]
    output_fname = args.output
    grid_fname   = args.grid_fname[0]
    
    # Read background & analysis
    #===========================
    bkg  = soca_cice5_utils.SeaIceState(fname=bkg_fname, gridname=grid_fname, ocnfname=ocn_fname, descriptor='Background')
    ana  = soca_cice5_utils.SeaIceState(fname=ana_fname, gridname=grid_fname, ocnfname=ocn_fname, descriptor='Analysis')
    incr = soca_cice5_utils.SeaIceState(gridname=grid_fname, descriptor='Increment')
    ana.info()

    #bkg.hist()
    #plt.show()

    #incr.delta(ana, bkg)
    #incr.update_ocn()
    incr.info()

    cmap=cm.jet
    #print 'Increment'
    #raw_input('<>?')
    #incr.plot(varname='aicen',vmin=-0.1,vmax=0.1,CMAP=cmap)
    #incr.plot(varname='hicen',vmin=-0.02,vmax=0.02,CMAP=cmap)
    #incr.plot(varname='hsnon',vmin=-0.01,vmax=0.01,CMAP=cmap)
    #incr.plot(varname='Tsfcn',vmin=-2.,vmax=2.,CMAP=cmap)
    #incr.plot(varname='',vmin=-0.2,vmax=0.2,CMAP=cmap)
    plt.show()    

    # Fix bounds for aicen, vicen and Tsfcn in analysis
    #==================================================
    ana.fix_bounds()

    # Define min/max & masking for incrememnt
    #========================================
    amin = 1e-3
    amax = 4.0
    tiny = 1e-3 # Masking parameter

    # Check aicen and fix bounds again
    #=================================
    for cat in range(bkg.ncat):   
        alfa=ana.aicen[cat,:,:]/bkg.aicen[cat,:,:]
        alfa[alfa>amax]  = amax
        alfa[alfa<amin]  = amin
        I=np.where(bkg.aicen[cat,:,:]<tiny)
        alfa[I]=1.0
        ana.aicen[cat,:,:]=alfa*bkg.aicen[cat,:,:]
    ana.fix_bounds(mask=True)
    bkg.info()
    ana.info()
    
    # Update ice and snow Volume, conserve background hicen and hsno
    #===============================================================
    for cat in range(bkg.ncat):    
        ana.vicen[cat,:,:]=ana.aicen[cat,:,:]*bkg.hicen[cat,:,:]
        ana.vsnon[cat,:,:]=ana.aicen[cat,:,:]*bkg.hsnon[cat,:,:]
    ana.fix_bounds()

    # Write analysis in clone of background file
    #===========================================
    shutil.copy(bkg_fname, output_fname)
    bkg.write(VAR=ana.aicen,fname=output_fname, varname='aicen')
    bkg.write(VAR=ana.vicen,fname=output_fname, varname='vicen')
    bkg.write(VAR=ana.vsnon,fname=output_fname, varname='vsnon')
    #bkg.write(VAR=ana.tsfcn,fname=output_fname, varname='Tsfcn')
    #bkg.write(VAR=ana.aicen,fname=output_fname, varname='aicen')

    # Save increment
    #===============
    incr.delta(ana, bkg)
    incr.info()

    incr.plot(varname='aicen',vmin=-0.025,vmax=0.025,figname='aicen.png')
    #incr.plot(varname='hicen',vmin=-0.02,vmax=0.02,figname='hicen.png')
    #incr.plot(varname='hsnon',vmin=-0.01,vmax=0.01,figname='hsnon.png')
    #incr.plot(varname='Tsfcn',vmin=-2.,vmax=2.,figname='Tsfcn.png')

    incr_fname='incr.nc'
    incr.write(VAR=incr.aicen,fname=incr_fname, varname='aicen',mode='w')
    incr.write(VAR=incr.vicen,fname=incr_fname, varname='vicen',mode='a')
    incr.write(VAR=incr.vsnon,fname=incr_fname, varname='vsnon',mode='a')
    incr.write(VAR=incr.tsfcn,fname=incr_fname, varname='Tsfcn',mode='a')
    incr.write(VAR=incr.aice,fname=incr_fname,  varname='aice',mode='a',size='2D')
    incr.write(VAR=incr.hice,fname=incr_fname,  varname='hice',mode='a',size='2D')
    incr.write(VAR=incr.dT,fname=incr_fname,  varname='dT',mode='a',size='2D')    

    #Save Ocean restart
    #T=soca_cice5_utils.readvar(ocn_fname,'Temp')
    #S=soca_cice5_utils.readvar(ocn_fname,'Salt')    
    #T[0,:,:]=T[0,:,:]+incr.dT
    
    #ncfile = netcdf.netcdf_file('ocn_rst_out.nc','w')
    #ncfile.createDimension('Time',1)
    #ncfile.createDimension('zaxis_1',40)
    #ncfile.createDimension('yaxis_1',410)
    #ncfile.createDimension('xaxis_1',720)    
    #VARtmp = ncfile.createVariable('Temp',np.dtype('float32').char,('Time','zaxis_1','yaxis_1','xaxis_1'))
    #VARtmp[0,:,:,:]=T
    #VARtmp = ncfile.createVariable('Salt',np.dtype('float32').char,('Time','zaxis_1','yaxis_1','xaxis_1'))
    #VARtmp[0,:,:,:]=S
    #ncfile.close()
    

    
    #I=np.where(bkg.aice<0.1)
    #AICE=reshape(ana.aicen[cat,I[0],I[1]],
    #plt.plot(ana.aicen[cat,I[0],I[1]]

    #plt.plot(np.squeeze(ana.aice.flatten()),np.squeeze(np.sum(ana.tsfcn,axis=0).flatten()),'.k')
    #plt.grid(True)
    #plt.show()

    #plt.figure(num=1)
    #plot2d(bkg.x, bkg.y, bkg.aicen, vmin=0.0, vmax=1.1)

    #plt.figure(num=2)
    #plot2d(ana.x, ana.y, ana.aicen, vmin=0.0, vmax=1.1)
    #plt.show()


