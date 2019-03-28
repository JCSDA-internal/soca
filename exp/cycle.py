#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Cycle SOCA

This script cycles mom6 through multiple DA windows.

Example:

     <---------------------- Window_length (24 hours) --------------------->

 ---|-----------|-----------|-----------|-----------|-----------|-----------|
    |                                   |                                   |      
 Window_begin: 2014-04-14T12:00:00Z     |                           2015-04-15T12:00:00Z
                                        |
                    datetime in obsfile: 2015-04-15T00:00:00Z


     $  ./cycle.py cycle.yml

Todo:

"""

import dateutil.parser
import datetime
from datetime import timedelta
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys, os
import yaml

def main(argv):
    description='Multiple DA Cycles. Ex: cycle.py cycle.yml'

    # Parse yaml config file
    configyml = sys.argv[1]
    with open(configyml, 'r') as f:
        config = yaml.load(f)

    ic_date    = str(config['Window_begin'])
    dt         = int(config['Window_length'])
    databases  = config['Observations']
    num_cycles = config['Num_cycles']
    
    ymdh = datetime.datetime.strptime(ic_date,"%Y%m%d%H")
    for t in range(num_cycles):
        # Extract year, month, day, hours from ic date
        yyyy = str(ymdh.year).zfill(4)
        mm   = str(ymdh.month).zfill(2)
        dd   = str(ymdh.day).zfill(2)
        hh   = str(ymdh.hour).zfill(2)

        # Get date for the middle of the window
        ymdh_mid = ymdh + timedelta(hours=int(dt))
        mid_date = ymdh_mid.strftime("%Y%m%d%H")

        # Loop through list of databases
        for database in databases:
            obs2ioda = database['Obs2ioda']
            binpath = database['binpath']
            databasepath = database['databasepath']
            
            command = binpath+obs2ioda+' --path '+databasepath+' --ic '+yyyy+mm+dd+hh+' --date '+mid_date+' --window '+str(dt)
            os.system(command)

        # Run DA
        command='./cycle.sh '+yyyy+mm+dd+hh+' '+str(dt)
        os.system(command)

        # Update time
        ymdh = ymdh + timedelta(hours=dt)

if __name__ == '__main__':
    main(sys.argv[1])
