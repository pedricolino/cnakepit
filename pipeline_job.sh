#!/bin/bash
#
#SBATCH --job-name=CNVsnake
#SBATCH --partition medium
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --export=all # Keep current environment variables
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cedric.moris@bih-charite.de
#SBATCH --output=logs/slurm_log/%x-%J.log

# Ensure logs folder exists -------------------------------------------------

mkdir -p logs/slurm_log/snakejobs
export SBATCH_DEFAULTS=" --output=logs/slurm_log/snakejobs/%x-%j.log"

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Activate bash cmd printing, debug info ------------------------------------

set -x
>&2 hostname
start_time=$(date +"%Y-%m-%d %H:%M:%S")
>&2 echo "Job started at $start_time"

# Activate conda environment ------------------------------------------------

# clutters log files and should test snakemake -n first anyway
#eval "$(conda shell.bash hook)"
#conda activate snakemake-vanilla

# Workflow specific parameters ----------------------------------------------

# For Mutect2 multithreading, see https://www.biostars.org/p/9549710/#9550707
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Kick off Snakemake --------------------------------------------------------

snakemake \
    --use-conda \
    --conda-frontend mamba \
    --cores 128 \
    --retries 2 \
    --profile cubi-v1 \
    --rerun-incomplete \
    --jobs=20 \
    --default-resources "runtime=3600" # 1 hour

# Finish up -----------------------------------------------------------------

end_time=$(date +"%Y-%m-%d %H:%M:%S")
>&2 echo "Job ended at $end_time"

# Calculate and print the total runtime in hours and minutes
start_seconds=$(date -d "$start_time" +"%s")
end_seconds=$(date -d "$end_time" +"%s")
runtime_seconds=$((end_seconds - start_seconds))

runtime_hours=$((runtime_seconds / 3600))
runtime_minutes=$(( (runtime_seconds % 3600) / 60 ))

>&2 echo "Total runtime: $runtime_hours hours and $runtime_minutes minutes"
>&2 echo "All done. Have a nice day."
