#!/bin/bash
#
#SBATCH --job-name=cnakepit
#SBATCH --partition medium
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --export=all # Keep current environment variables
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cedric.moris@bih-charite.de
#SBATCH --output=logs/slurm_log/%x-%J.log

start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Starting CNAkepit pipeline at $start_time. Look out for snakes..."

# Ensure logs folder exists -------------------------------------------------

mkdir -p logs/slurm_log/snakejobs
export SBATCH_DEFAULTS=" --job-name {rule}.{wildcards} --output=logs/slurm_log/snakejobs/%x-%j.log"

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Activate bash cmd printing, debug info ------------------------------------

set -x
>&2 hostname

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
    --retries 2 \
    --profile cubi-v1 \
    --rerun-incomplete \
    --jobs=200 \
    --slurm \
    --default-resources slurm_account=hpc-ag-cubi slurm_partition=short "runtime=240"
#    --batch cnvkit_heatmap_hmm=1/6 # does not work well bc rule autobin requires all files/samples, quite early on in the pipeline


# Finish up -----------------------------------------------------------------

set +x

end_time=$(date +"%Y-%m-%d %H:%M:%S")
>&2 echo "Job ended at $end_time"

# Calculate and print the total runtime in hours and minutes
start_seconds=$(date -d "$start_time" +"%s")
end_seconds=$(date -d "$end_time" +"%s")
runtime_seconds=$((end_seconds - start_seconds))

runtime_hours=$((runtime_seconds / 3600))
runtime_minutes=$(( (runtime_seconds % 3600) / 60 ))

>&2 echo "Total runtime: $runtime_hours hours and $runtime_minutes minutes"
>&2 echo "All done. You may have leave the cnakepit now. Have a nice day."
