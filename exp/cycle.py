#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Cycle SOCA

This scripts cycles mom6 through multiple DA windows.

Example:
     $  ./_.................

Todo:

"""

import dateutil.parser
import datetime
from datetime import timedelta
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys, os

if __name__ == '__main__':

    description='Multiple DA Cycles. Ex: cycle.py --ic 2018071912 --window 24 --num_cycles 10'
    parser = ArgumentParser(description=description,formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ic', help='initial condition, yyyymmddhh', type=str, required=True)    
    parser.add_argument('--window', help='half window length in hours', type=str, required=True)
    parser.add_argument('--num_cycles', help='number of DA cycles', type=str, required=True)    
    args = parser.parse_args()

    dt = int(args.window)
    ic_date = args.ic
    num_cycles = int(args.num_cycles)

    ymdh = datetime.datetime.strptime(ic_date,"%Y%m%d%H")
    for t in range(num_cycles):    
        yyyy=str(ymdh.year).zfill(4)
        mm=str(ymdh.month).zfill(2)
        dd=str(ymdh.day).zfill(2)
        hh=str(ymdh.hour).zfill(2)
        command='cycle.sh '+yyyy+' '+mm+' '+dd+' '+hh+' '+str(dt)
        print command
        os.system(command)
        ymdh = ymdh + timedelta(hours=dt)
