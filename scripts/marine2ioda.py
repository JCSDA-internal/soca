#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import yaml
from glob import glob
from datetime import datetime, timedelta
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Import marine converters
import godae_profile2ioda
import godae_trak2ioda
import godae_ship2ioda
import rads_adt2ioda
import smap_sss2ioda
import gds2_sst2ioda


def create_datedict(adate, twin=24):

    datedict = {}

    dates_ymd = []
    cdate = adate - twin/2
    while cdate <= adate + twin/2:
        dates_ymd.append(cdate.strftime('%Y%m%d'))
        cdate += timedelta(days=1)

    datedict['analysis_date'] = adate
    datedict['analysis_window_hrs'] = twin
    datedict['obs_dates_ymd'] = dates_ymd

    return datedict


def godae(obs, obspaths, datedict):

    otype = obs['obstype']

    print('processing ... %s' % otype)

    datain = obs.get('datain', obspaths.get('datain'))
    dataout = obs.get('dataout', obspaths.get('dataout'))

    if otype in ['godae_prof']:
        fext = 'profile'
    elif otype in ['godae_ship']:
        fext = 'ship'
    elif otype in ['godae_trak']:
        fext = 'trak'
    else:
        raise KeyError('Unknown observation type %s' % otype)

    files = []
    dates_ymd = datedict['obs_dates_ymd']
    for cdate in dates_ymd:
        datapath = '%s/%s/%s/*.%s' % (datain, obs['source'], cdate, fext)
        files.append(glob(datapath))

    flatten = lambda l: [item for sublist in l for item in sublist]
    files = flatten(files)

    if len(files) == 0:
        print('No files to process for %s, RETURN!' % otype)
        return

    adatestr = datedict['analysis_date'].strftime('%Y%m%d%H')
    fout = '%s/%s.%s.nc' % (dataout, adatestr, otype)
    try:
        os.makedirs(os.path.dirname(fout))
    except:
        print('%s exists, use it!' % os.path.dirname(fout))

    sys.argv = ['%s' % obs['converter'],
                '-d', '%s' % adate.strftime('%Y%m%d%H'),
                '-o', '%s' % (fout),
                '-i']
    sys.argv += files

    print('running: %s' % ' '.join(sys.argv))

    if otype in ['godae_prof']:
        getattr(godae_profile2ioda, 'main')()
    elif otype in ['godae_ship']:
        getattr(godae_ship2ioda, 'main')()
    elif otype in ['godae_trak']:
        getattr(godae_trak2ioda, 'main')()

    return


def adt(obs, obspaths, datedict):

    otype = obs['obstype']

    print('processing ... %s' % otype)

    datain = obs.get('datain', obspaths.get('datain'))
    dataout = obs.get('dataout', obspaths.get('dataout'))

    files = []
    dates_ymd = datedict['obs_dates_ymd']
    for cdate in dates_ymd:
        datapath = '%s/%s/%s/*' % (datain, obs['source'], cdate)
        files.append(glob(datapath))

    flatten = lambda l: [item for sublist in l for item in sublist]
    files = flatten(files)

    if len(files) == 0:
        print('No files to process for %s, RETURN!' % otype)
        return

    adatestr = datedict['analysis_date'].strftime('%Y%m%d%H')
    fout = '%s/%s.%s.nc' % (dataout, adatestr, otype)
    try:
        os.makedirs(os.path.dirname(fout))
    except:
        print('%s exists, use it!' % os.path.dirname(fout))

    sys.argv = ['%s' % obs['converter'],
                '-d', '%s' % adate.strftime('%Y%m%d%H'),
                '-o', '%s' % (fout),
                '-i']
    sys.argv += files

    print('running: %s' % ' '.join(sys.argv))

    getattr(rads_adt2ioda, 'main')()

    return


def sss(obs, obspaths, datedict):

    otype = obs['obstype']

    print('processing ... %s' % otype)

    datain = obs.get('datain', obspaths.get('datain'))
    dataout = obs.get('dataout', obspaths.get('dataout'))

    files = []
    dates_ymd = datedict['obs_dates_ymd']
    for cdate in dates_ymd:
        datapath = '%s/%s/%s/*' % (datain, obs['source'], cdate)
        files.append(glob(datapath))

    flatten = lambda l: [item for sublist in l for item in sublist]
    files = flatten(files)

    if len(files) == 0:
        print('No files to process for %s, RETURN!' % otype)
        return

    adatestr = datedict['analysis_date'].strftime('%Y%m%d%H')
    fout = '%s/%s.%s.nc' % (dataout, adatestr, otype)
    try:
        os.makedirs(os.path.dirname(fout))
    except:
        print('%s exists, use it!' % os.path.dirname(fout))

    sys.argv = ['%s' % obs['converter'],
                '-d', '%s' % adate.strftime('%Y%m%d%H'),
                '-o', '%s' % (fout),
                '-i']
    sys.argv += files

    print('running: %s' % ' '.join(sys.argv))

    getattr(smap_sss2ioda, 'main')()

    return


def sst(obs, obspaths, datedict):

    otype = obs['obstype']

    print('processing ... %s' % otype)

    datain = obs.get('datain', obspaths.get('datain'))
    dataout = obs.get('dataout', obspaths.get('dataout'))

    platform = obs.get('platform', ['*'])
    threads = obs.get('threads', 1)
    thin = obs.get('thin', 0.)

    files = []
    dates_ymd = datedict['obs_dates_ymd']
    for cdate in dates_ymd:
        for plat in platform:
            datapath = '%s/%s/%s/*%s*' % (datain, obs['source'], cdate, plat)
            files.append(glob(datapath))

    flatten = lambda l: [item for sublist in l for item in sublist]
    files = flatten(files)

    if len(files) == 0:
        print('No files to process for %s, RETURN!' % otype)
        return

    adatestr = datedict['analysis_date'].strftime('%Y%m%d%H')
    fout = '%s/%s.%s.nc' % (dataout, adatestr, otype)
    try:
        os.makedirs(os.path.dirname(fout))
    except:
        print('%s exists, use it!' % os.path.dirname(fout))

    sys.argv = ['%s' % obs['converter'],
                '--threads', '%d' % threads,
                '--thin', '%f' % thin,
                '-d', '%s' % adate.strftime('%Y%m%d%H'),
                '-o', '%s' % (fout),
                '-i']
    sys.argv += files

    print('running: %s' % ' '.join(sys.argv))

    getattr(gds2_sst2ioda, 'main')()

    return


if __name__ == '__main__':

    desc = 'Convert raw marine observations into IODA netCDF4 format'
    parser = ArgumentParser(
        description=desc,
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--input', help='name of the YAML configuration file',
        type=str, required=True)
    parser.add_argument(
        '-d', '--date', help='analysis date',
        type=str, metavar='YYYYMMDDHH', required=True)
    parser.add_argument(
        '-w', '--window', help='analysis window (hrs)',
        type=int, default=24, required=False)

    args = parser.parse_args()
    ymlfile = args.input
    adate = datetime.strptime(args.date, '%Y%m%d%H')
    twin = timedelta(hours=args.window)

    with open(ymlfile, 'r') as fh:
        ymlcfg = yaml.safe_load(fh)

    datedict = create_datedict(adate, twin)

    obspaths = ymlcfg['obspaths']
    obs2ioda = ymlcfg['obs2ioda']

    for o in obs2ioda:

        otype = o['obstype']

        if otype.startswith('godae'):
            godae(o, obspaths, datedict)

        elif otype.startswith('sst'):
            sst(o, obspaths, datedict)

        elif otype in ['adt']:
            adt(o, obspaths, datedict)

        elif otype in ['sss']:
            sss(o, obspaths, datedict)

    sys.exit(0)
