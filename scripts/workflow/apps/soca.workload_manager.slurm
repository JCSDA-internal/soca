function submitjob {
    sbatch --ntasks $JOB_NPES --time=$JOB_TIME --qos=debug -A $JOB_ACCT -J $JOB_NAME \
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
