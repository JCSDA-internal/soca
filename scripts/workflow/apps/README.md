The `cycle.sh` script in this directory can be used to run cycling experiments on any HPC
system that uses SLURM or singularity containers. It implements a very
simple strategy of using a single job to
run all the steps of the cycle. At the end of a cycle the script checks to see if there is
is enough time remaining to keep going, and if there is not enough time, the script resubmits itself
to the job queue and terminates.

### Preparing ###

The following assumptions are made (you're on your own for these parts):
* MOM6/SIS2 standalone executable has been compiled separately
* SOCA has been compiled
* binary files required by various SOCA/MOM configurations have been downloaded from the
shared google drive [mom6sis2_global.tgz](https://drive.google.com/open?id=1uDwWsuaMs_8M8yHiIVJRhGj99dYtDSKX)
* preprocessed CFSR surface forcing files are already present
* observations have already been converted into ioda format

### Setting up a new experiment ###

Create a directory for the experiment and copy the slurm cycling files

```
mkdir expdir
cd expdir
cp <soca_dir>/scripts/workflow/apps/* .
```

Edit the experiment configuration file `exp.config` paying close attention to make sure the
default direcotry locations have been changed to something appropriate for you.

run `cycle.h` from the experiment directory. This will launch a slurm job, and if everything goes
well, will run the desired duration of the experiment and create the following output directories (note: all netcdf files
will be compressed. State-space files also have their precision reduced before being compressed)
* `ana/` -  analysis files
* `bkg/` -  background files
* `da_init/` - files required by the DA (bump/geometry) that are created on the first cycle, and used as is after that
* `diag/` - diagnostic files from the MOM6 forecast
* `logs/` - logs, the overall cycle output is placed in `log/slurm.log`. All individual job steps are placed in their own log file
* `omb/` - observation space statistics from the DA output
* `rst/` - MOM restart file that are saved at the end of a cycle for use at the beginning of
  the next cycle. Only the most recent cycle is saved
