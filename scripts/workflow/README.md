The scripts in this directory can be used for cycling experiments of SOCA. For directions on running check
the directory of the workflow system you want to use (only option right now is SLURM or none).

* `subscripts/` - contains the scripts for the individual job steps that should hopefully be agnostic
of the system being run on and the workflow type being used

## workflow system ##
* `apps/` - cycling script that should work on any HPC that uses the slurm job
managemnt system or workstation with no workload management.
