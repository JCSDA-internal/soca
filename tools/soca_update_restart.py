#!/usr/bin/env python
"""
SOCA Restart File Updater

This script updates SOCA model restart files with analysis data from specified netCDF files, 
as directed by a YAML configuration file. It applies specified bounds checks and maximum increment limits 
to the analysis data before updating the restart files.

Usage:
    ./soca_restart_updater.py <input_yaml_configuration_file>

Input Configuration File Format:
    The configuration file should be in YAML format, specifying input and output restart files, 
    and detailing the analysis files along with variable names and constraints. An example configuration:

    restarts:
    - input restart: path/to/input/MOM.res.nc
      output restart: path/to/output/test.nc
      input analysis:
      - file: path/to/analysis/ocn.enspert.ens.1.2018-04-15T00:00:00Z.PT0S.nc
        variables:
        - name: Temp
          min: -1.8
          max: 32.0
        - name: u
          increment max: 0.5
        [...]

Features:
    - Copies input restart files and updates them with analysis data.
    - Applies constraints on the analysis data, such as minimum and maximum values, and maximum increments.
    - Validates configuration file for required fields and logical consistency of specified constraints.
"""

import yaml
import shutil
import netCDF4 as nc
import numpy as np
import argparse
import os

# ---------------------------------------------------------------------------------------
def validate_config(config):
  # error checking, curtesy of chatgpt
  if 'restarts' not in config or not config['restarts']:
    raise ValueError("Configuration must contain at least one 'restarts' item")
  for restartConfig in config['restarts']:
    if 'input restart' not in restartConfig or 'output restart' not in restartConfig or 'input analysis' not in restartConfig:
      raise ValueError("Missing required keys in restart configuration")
    if not restartConfig['input analysis']:
      raise ValueError("'input analysis' must contain at least one item")
    for anaConfig in restartConfig['input analysis']:
      if 'file' not in anaConfig or 'variables' not in anaConfig:
        raise ValueError("Missing 'file' or 'variables' keys in input analysis configuration")
      if not anaConfig['variables']:
        raise ValueError("'variables' must contain at least one item")
      for varConfig in anaConfig['variables']:
        if 'name' not in varConfig:
          raise ValueError("Missing 'name' key in variable configuration")
        if 'min' in varConfig and 'max' in varConfig:
          if varConfig['min'] >= varConfig['max']:
            raise ValueError("'min' value must be less than 'max' value for variable")

def file_exists(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
# ---------------------------------------------------------------------------------------

# get input configuration
parser = argparse.ArgumentParser(description="SOCA utility to update a full restart file with an analysis file.")
parser.add_argument("INPUT", help="input yaml configuration file")
args = parser.parse_args()
with open(args.INPUT, 'r') as config_file:
  config = yaml.safe_load(config_file)
validate_config(config)

# we're going to modify 1 or more restart files
for restartConfig in config['restarts']:

  # make a copy of the input file
  inputRestart = restartConfig['input restart']
  outputRestart =  restartConfig['output restart'] + '.tmp'
  file_exists(inputRestart)
  shutil.copy(inputRestart, outputRestart)
  print ("Creating ", restartConfig['output restart'])

  # open the files and process each input analysis 
  with nc.Dataset(outputRestart, 'r+') as restartNcd:
    for anaConfig in restartConfig['input analysis']:
      file_exists(anaConfig['file'])
      with nc.Dataset(anaConfig['file'], 'r') as inputNcd:
        
        # for each input variable
        for varConfig in anaConfig['variables']:
          var_name = varConfig['name']
          print("  variable: ", var_name)

          ana_var = inputNcd.variables[var_name][:]
          rst_var = restartNcd.variables[var_name][:]

          # bounds check on increment
          if 'increment max' in varConfig:
            maxInc = varConfig['increment max']
            print("    imposing max increment of ", maxInc)
            inc = ana_var - rst_var
            inc = np.clip(inc, -maxInc, maxInc)
            ana_var = rst_var + inc

          # bounds check on analysis
          if 'min' in varConfig:
            val = varConfig['min']
            print("    imposing min value of ", val)
            ana_var[ana_var < val] = val
          if 'max' in varConfig:
            val = varConfig['max']
            print("    imposing max value of ", val)
            ana_var[ana_var > val] = val

          # write out
          restartNcd.variables[var_name][:] = ana_var

  # all done with the restart file, move from .tmp to final filename
  shutil.move(outputRestart, restartConfig['output restart'])
  print ("")