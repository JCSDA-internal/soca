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
    bkg  = soca_cice5_utils.SeaIceState(fname=bkg_fname,
                                        gridname=grid_fname,
                                        ocnfname=ocn_fname,
                                        descriptor='Background')
    ana  = soca_cice5_utils.SeaIceState(fname=ana_fname,
                                        gridname=grid_fname,
                                        ocnfname=ocn_fname,
                                        descriptor='Analysis')
    bkg.info()
    ana.info()

    # Fix bounds for aicen, vicen and Tsfcn in analysis
    #==================================================
    ana.fix_bounds()
    ana.update()
    ana.info()

    # Define min/max & masking for incrememnt
    #========================================
    amin = 1e-6
    amax = 10.0
    tiny = 1e-6 # Masking parameter

    alfa = np.ones(np.shape(bkg.aice))
    I=np.where(bkg.aice>tiny)
    alfa[I]=ana.aice[I]/bkg.aice[I]
    alfa[alfa>amax]  = amax
    alfa[alfa<amin]  = amin

    alfa=0.9999999*alfa
    # Update ice and snow
    #====================
    ana.aicen=alfa*bkg.aicen
    ana.vicen=alfa*bkg.vicen
    ana.vsnon=alfa*bkg.vsnon
    ana.update()
    ana.info()

    # Write analysis in clone of background file
    #===========================================
    shutil.copy(bkg_fname, output_fname)
    bkg.write(VAR=ana.aicen,fname=output_fname, varname='aicen')
    bkg.write(VAR=ana.vicen,fname=output_fname, varname='vicen')
    bkg.write(VAR=ana.vsnon,fname=output_fname, varname='vsnon')


