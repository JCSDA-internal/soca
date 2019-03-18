#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Cycle SOCA

This script cycles mom6 through multiple DA windows.

Example:

     <---------------------- window_length (2 hours) ---------------------->

 ---|-----------|-----------|-----------|-----------|-----------|-----------|
    |                                   |                                   |      
 window_begin: 2014-04-14T12:00:00Z     |                           2015-04-15T12:00:00Z
                                        |
                    datetime in obsfile: 2015-04-15T00:00:00Z


     $  ./cycle.py --ic 2018041412 --window 24 --num_cycles 10

Todo:

"""

import dateutil.parser
import datetime
from datetime import timedelta
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys, os

if __name__ == '__main__':
    description='Multiple DA Cycles. Ex: cycle.py --ic 2018041500 --window 24 --num_cycles 2'
    parser = ArgumentParser(description=description,formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ic', help='initial condition, yyyymmddhh', type=str, required=True)    
    parser.add_argument('--window', help='window length in hours', type=str, required=True)
    parser.add_argument('--num_cycles', help='number of DA cycles', type=str, required=True)    
    args = parser.parse_args()

    dt = int(args.window)
    ic_date = args.ic
    num_cycles = int(args.num_cycles)

    ymdh = datetime.datetime.strptime(ic_date,"%Y%m%d%H")
    for t in range(num_cycles):
        # Extract year, month, day, hours from date
        yyyy=str(ymdh.year).zfill(4)
        mm=str(ymdh.month).zfill(2)
        dd=str(ymdh.day).zfill(2)
        hh=str(ymdh.hour).zfill(2)

        # Get date for the middle of the window
        ymdh_mid = ymdh + timedelta(hours=dt)
        
        # Prepare ADT observation files
        adt_script = '../../bin/altimeter2ioda.py'
        adt_path = '../../../../../../Data/nesdis/rads/'
        mid_date = str(ymdh_mid.year).zfill(4)+str(ymdh_mid.month).zfill(2)+str(ymdh_mid.day).zfill(2)+str(ymdh_mid.hour).zfill(2)
        command=adt_script+' --path '+adt_path+' --ic '+yyyy+mm+dd+hh+' --date '+mid_date+' --window '+str(dt)
        print command
        os.system(command)
        
        # Run DA
        command='cycle.sh '+yyyy+' '+mm+' '+dd+' '+hh+' '+str(dt)
        print command
        os.system(command)

        # Update time
        ymdh = ymdh + timedelta(hours=dt)
