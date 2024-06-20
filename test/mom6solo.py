#!/usr/bin/env python3
import sys
import os
import yaml
import subprocess
import shutil
import re
import pandas as pd
from datetime import datetime, timedelta

MOM6EXE = 'mom6solo'

# get input parameters
binDir = os.environ['BIN_DIR']
mpiExe = os.environ['MPIEXE']
with open(sys.argv[1], 'r') as config_file:
  config = yaml.safe_load(config_file)

#----------------------------------------------------------------------------------------
# prepare input.nml
#----------------------------------------------------------------------------------------
fcst_len = round(pd.Timedelta(config['forecast length']).total_seconds() / 3600)
replacements = {
  'ocean_solo_nml' : {
    'months' : '0',
    'days' : '0',
    'hours': fcst_len
  },
}

# read input file
with open(config['mom6_input_nml'], 'r') as inputFile:
  inputNml = inputFile.read()

# do the replacements
def replace_config_content(content, replacements):
  for section, changes in replacements.items():
    section_pattern = rf"(&{section}\s*.*?)/\s*\n"
    section_match = re.search(section_pattern, content, re.DOTALL)
    if section_match:
      section_content = section_match.group(0)
      for key, value in changes.items():
        def replacement_function(match):
          start, end = match.group(1), match.group(3)
          return f"{start}{value}{end}"                
        key_pattern = rf"(\b{key}\s*=\s*)([^,\n]*)(,|/|\n)"
        section_content = re.sub(key_pattern, replacement_function, section_content, flags=re.MULTILINE)
      # Replace the old section content with the new one
      content = content.replace(section_match.group(0), section_content)
  return content
modified_inputNml = replace_config_content(inputNml, replacements)

# write output file
with open('input.nml', 'w') as file:
    file.write(modified_inputNml)

#----------------------------------------------------------------------------------------
# prepare MOM_override (for setting intermediate restart file frequency)
#----------------------------------------------------------------------------------------
restartFrequency = round(pd.Timedelta(config['output']['frequency']).total_seconds() / 3600)
overrides = {
  'RESTART_CONTROL' : 2,
  'RESTINT' : restartFrequency / 24.0
}
with open('MOM_override', 'w') as file:
  for key, value in overrides.items():
    file.write(f"#override {key} = {value}\n")


#----------------------------------------------------------------------------------------
# link input restart files
#----------------------------------------------------------------------------------------
if os.path.exists ('RESTART'):
  shutil.rmtree('RESTART')  
os.makedirs('RESTART')
if os.path.exists ('RESTART_IN'):
  shutil.rmtree('RESTART_IN')  
os.makedirs('RESTART_IN')
os.symlink("../"+config['initial condition']['ocn_filename'], 'RESTART_IN/MOM.res.nc')


#----------------------------------------------------------------------------------------
# prepare ocean_solo.res (for setting initial date)
#----------------------------------------------------------------------------------------
init_date = config['initial condition']['date']
with open('RESTART_IN/ocean_solo.res', 'w') as file:
  file.write("     4        (Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)\n")
  file.write(f"  {init_date.year:4d}     {init_date.month:2d}    {init_date.day:2d}     {init_date.hour:2d}     {init_date.minute:2d}     {init_date.second:2d}        Model start time:   year, month, day, hour, minute, second\n")
  file.write(f"  {init_date.year:4d}     {init_date.month:2d}    {init_date.day:2d}     {init_date.hour:2d}     {init_date.minute:2d}     {init_date.second:2d}        Current model time: year, month, day, hour, minute, second\n")
 

#----------------------------------------------------------------------------------------
# run executable
#----------------------------------------------------------------------------------------
exePath = binDir + '/' + MOM6EXE
try:
  # Capture standard output and error
  result = subprocess.run(mpiExe+" "+exePath, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
  # If the command was successful, print the output
  print(f"Executable ran successfully: {result.stdout}")
except subprocess.CalledProcessError as e:
  # If the command failed, print the error output and raise an exception
  print(f"Error running executable: {e.stderr}")
  raise Exception(f"Executable failed with return code {e.returncode}")


#----------------------------------------------------------------------------------------
# Move/rename the output restart files
#----------------------------------------------------------------------------------------
for fcstHr in range(restartFrequency, fcst_len+1, restartFrequency):
  rstDate = init_date + timedelta(hours=fcstHr)
  
  # get rstDate in format such as MOM.res_Y2018_D105_S03600.nc
  rst_day = rstDate.timetuple().tm_yday  
  rst_seconds = rstDate.hour * 3600
  rst_src = f"MOM.res_Y{rstDate.year}_D{rst_day:03}_S{rst_seconds:05}.nc"
 
  # get init_date / fcstHr in a format such as ocn.forecast_mom6.fc.2018-04-15T00:00:00Z.PT1H.nc
  forecast_date = init_date + timedelta(hours=fcstHr)
  forecast_date_str = init_date.strftime("%Y-%m-%dT%H:%M:%SZ")
  rst_dst = f"ocn.{config['output']['exp']}.fc.{forecast_date_str}.PT{fcstHr}H.nc"
    
  shutil.move('RESTART/'+rst_src, config['output']['datadir']+'/'+rst_dst)
