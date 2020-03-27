#!/bin/bash
#================================================================================
# SLURM cycling script for SOCA experiments.
#
# A very simple workflow that calls all steps of the DA cycle sequentially.
# Run this script directly from the experiment directory to start the job submission.
#
# This script will resubmit itself to SLURM and exit if there is not enough
# time remaining for the current job (based on an averge cycle runtime kept track
# of by this script)
#================================================================================
set -eu

# read the experiment configuration file
source exp.config

# make sure logging directory is setup
LOG_DIR=$EXP_DIR/logs
mkdir -p $LOG_DIR

# are we under SLURM? If not resubmit this script under sbatch
function submitjob {
    sbatch -N $JOB_NODES --time=$JOB_TIME -A $JOB_ACCT -J $JOB_NAME \
	   -o $LOG_DIR/slurm.log \
	   $EXP_DIR/cycle.sh    
}
set +u
if [[ -z "$SLURM_JOB_ID" ]]; then
    set -u
    echo "launching under SLURM..."
    submitjob
    exit 0
fi
set -u


#================================================================================
# Start of loop
#================================================================================
export TIMEFORMAT='%1Rs'
cycle_avg_count=0
cycle_avg_runtime=0
while true; do
    cycle_start=$(date +%s)
    
    # determine where the experiment left off and if this is the first cycle
    CYCLE_STATUS_FILE=$EXP_DIR/rst/cycle_status
    if [[ ! -e $CYCLE_STATUS_FILE ]]; then
	FCST_START_TIME=$(date -ud "$EXP_START_TIME")
	FCST_RESTART=0
    else
	FCST_RESTART=1
	FCST_START_TIME=$(cat $CYCLE_STATUS_FILE)
    fi

    # are we done with the experiment?
    if [[ $(date -ud "$FCST_START_TIME" +%Y%m%d%H) -ge \
	  $(date -ud "$EXP_END_TIME" +%Y%m%d%H) ]]; then
	echo "Done with the experiment."
	exit 0
    fi
    echo -e "\n\nStarting cycle for $FCST_START_TIME"

    # resubmit this job if we are almost out of time.
    # (don't check this until cycle has run at least once, obviously)
    if [[ $cycle_avg_count -gt 0 ]]; then
       end_time=$(squeue -h -j $SLURM_JOB_ID -o %e)
       end_time=$(date -d "${end_time:0:10} ${end_time:11:8}" +%s)
       ((rem_time=end_time-cycle_start))
       min_time=$(bc <<< "$cycle_avg_runtime * 1.2 / 1") # fudge factor of 1.2x for making sure
                                                         # there is enough runtime left
       echo "$rem_time seconds remaining for current SLURM job"
       if [[ "$rem_time" -lt "$min_time" ]]; then
	   echo "Almost out of time, resubmitting a new job"
	   submitjob
	   exit 0
       fi
    fi
    
    # determine other derived variables    
    FCST_END_TIME=$(date -ud "$FCST_START_TIME + $FCST_LEN hours")
    ANA_TIME="$FCST_END_TIME"    
    SCRATCH_DIR_CYCLE=$SCRATCH_DIR/$(date -ud "$ANA_TIME" +%Y%m%d%H)
    LOG_DIR_CYCLE=$LOG_DIR/$(date -ud "$ANA_TIME" +%Y%m%d%H)
    SCRIPT_DIR=$SOCA_BUNDLE_DIR/soca/scripts/workflow/subscripts
    
    # setup working/output directories
    mkdir -p $EXP_DIR/rst
    mkdir -p $SCRATCH_DIR_CYCLE
    mkdir -p $LOG_DIR_CYCLE
    

    #------------------------------------------------------------
    # fcst.prep.sh
    #------------------------------------------------------------
    (
	echo -n "Running fcst.prep.sh ... "
	export ENS_SIZE=1
	export FCST_END_TIME=$(date -ud "$FCST_END_TIME" +%Y%m%d)
	export FCST_START_TIME=$(date -ud "$FCST_START_TIME" +%Y%m%d)
	export FORC_MEAN_FILE
	export FORC_RUNOFF=0
	export FORC_VAR
	export FORC_VAR_POS=" "
	export IC_GEN=0
	export WORK_DIR=$SCRATCH_DIR_CYCLE/fcst.prep
	time $SCRIPT_DIR/fcst.prep.sh &> $LOG_DIR_CYCLE/fcst.prep
    ) || { echo "ERROR in fcst.prep"; exit 1; }

    
    #------------------------------------------------------------
    # fcst.run.sh
    #------------------------------------------------------------
    (
	echo -n "Running fcst.run.sh ... "	
	export FCST_LEN
	export FCST_RESTART
	export FCST_START_TIME
	export FORC_ATM
	export FORC_DIR=$SCRATCH_DIR_CYCLE/fcst.prep/forc/mem_0000
	export MOM_CONFIG=$MODEL_CONFIG/model
	export MOM_DATA=$MODEL_DATA/model
	export MOM_EXE
	export MOM_IC
	export RESTART_DIR_IN=$EXP_DIR/rst/$(date -ud "$FCST_START_TIME" +%Y%m%d%H)
	export WORK_DIR=$SCRATCH_DIR_CYCLE/fcst.run
	time $SCRIPT_DIR/fcst.run.sh &> $LOG_DIR_CYCLE/fcst.run
    ) || { echo "ERROR in fcst.run"; exit 1; }    

    #------------------------------------------------------------
    # da.init.sh
    #------------------------------------------------------------
    (
	echo -n "Running da.init.sh ... "
	export DA_INIT_DIR=$EXP_DIR/da_init
	export FCST_LEN
	export FCST_START_TIME
	export MOM_CONFIG=$MODEL_CONFIG/model
	export MOM_DATA=$MODEL_DATA/model
	export RESTART_DIR=$SCRATCH_DIR_CYCLE/fcst.run/RESTART
	export SOCA_BIN_DIR
	export SOCA_CONFIG=$MODEL_CONFIG/soca
	export SOCA_DATA=$MODEL_DATA/soca
	export WORK_DIR=$SCRATCH_DIR_CYCLE/da.init
	time $SCRIPT_DIR/da.init.sh &> $LOG_DIR_CYCLE/da.init
    ) || { echo "ERROR in da.init"; exit 1; }

    #------------------------------------------------------------
    # da.run.sh
    #------------------------------------------------------------
    (
	echo -n "Running da.run.sh ... "
	export ANA_TIME
	export CYCLE_RST_DIR=$EXP_DIR/rst/$(date -ud "$ANA_TIME" +%Y%m%d%H)
	export DA_INIT_DIR=$EXP_DIR/da_init
	export FCST_LEN
	export FCST_START_TIME
	export MOM_CONFIG=$MODEL_CONFIG/model
	export MOM_DATA=$MODEL_DATA/model	
	export OBS_IODA
	export RESTART_DIR=$SCRATCH_DIR_CYCLE/fcst.run/RESTART
	export SOCA_BIN_DIR
	export SOCA_CONFIG=$MODEL_CONFIG/soca
	export SOCA_DATA=$MODEL_DATA/soca
	export WORK_DIR=$SCRATCH_DIR_CYCLE/da.run
	time $SCRIPT_DIR/da.run.sh &> $LOG_DIR_CYCLE/da.run
    ) || { echo "ERROR in da.run"; exit 1; }

    
    #------------------------------------------------------------
    # da.post.sh
    #------------------------------------------------------------
    (
	echo -n "Running da.post.sh ... "
	export ANA_TIME
	export OUT_DIR=$EXP_DIR
	export DA_ANA=$SCRATCH_DIR_CYCLE/da.run/RESTART/MOM.res.nc
	export DA_BKG=$SCRATCH_DIR_CYCLE/da.run/INPUT/MOM.res.nc
	export DA_OMB=$SCRATCH_DIR_CYCLE/da.run/obs_out
	export FCST_DIAG=$SCRATCH_DIR_CYCLE/fcst.run/*.ocean_diag.nc
	export WORK_DIR=$SCRATCH_DIR_CYCLE/da.post
	time $SCRIPT_DIR/da.post.sh &> $LOG_DIR_CYCLE/da.post
    ) || { echo "ERROR in da.post"; exit 1; }

    
    #============================================================
    # done with this day of the cycle, cleanup and prepare for the next cycle
    echo "$ANA_TIME" > $CYCLE_STATUS_FILE
    rm -r $SCRATCH_DIR_CYCLE
    if [[ $FCST_RESTART == 1 ]]; then
	rm -r $EXP_DIR/rst/$(date -ud "$FCST_START_TIME" +%Y%m%d%H)
    fi

    # update statistics on average cycle runtime
    cycle_end=$(date +%s)
    ((cycle_runtime=cycle_end-cycle_start))
    ((cycle_avg_runtime=(cycle_avg_runtime*cycle_avg_count+cycle_runtime)/(cycle_avg_count+1) ))
    ((cycle_avg_count=cycle_avg_count+1))
    echo "Cycle runtime:     $cycle_runtime seconds"
    echo "Cycle runtime avg: $cycle_avg_runtime seconds"

done
