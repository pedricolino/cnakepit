#!/bin/bash
#
#SBATCH --job-name=MH-1000-cnakepit
#SBATCH --partition medium
#SBATCH --mem=32G
#SBATCH --ntasks=8
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --export=all # Keep current environment variables
#SBATCH --output=logs/slurm_log/%x-%J.log
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=your@email.address

start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Starting CNAkepit pipeline at $start_time. Look out for snakes..."

#--- Activate bash cmd printing, debug info ------------------------------------

set -x
hostname >&2

#--- Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

#--- Activate conda environment ------------------------------------------------

## Possible but clutters log files & should test w/ snakemake -n first anyway
# eval "$(conda shell.bash hook)"
# conda activate snakemake-vanilla

#--- Workflow specific parameters ----------------------------------------------

# For Mutect2 multithreading, see https://www.biostars.org/p/9549710/#9550707
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#--- Kick off Snakemake --------------------------------------------------------

# Run the pipeline for 10 different cohortsizes, each 10 times
for run_number in {1..100}; do
    for cohort_size in 50 100 300 600; do
        snakemake \
            --use-conda \
            --retries 1 \
            --conda-frontend conda \
            --conda-prefix ~/work/miniconda/envs \
            --profile cubi-v1 \
            --rerun-incomplete \
            --rerun-triggers mtime \
            --jobs 150 \
            --cores 1800 \
            --slurm \
            --max-jobs-per-second 1 \
            --keep-going \
            --default-resources slurm_account=hpc-ag-cubi slurm_partition=short "runtime=240" \
            --config "pon={\"reference_sample_names\":\"results/stats/low_ploidy_PONs_for_different_cohortsizes_connected_runs/lowploidy_samples_${cohort_size}_run${run_number}.tsv\", \"directory_suffix\":\"lowploidy_samples_${cohort_size}_run${run_number}_connected_runs\"}"
    done
done

#    --batch cnvkit_heatmap_hmm=1/6 # process only 1/6 of the samples at a time, up to the rule that requires all samples
#   --rerun-triggers mtime # prevents rerunning of rules because of (minor) code changes
#   --keep-going # continue running even if one or few samples continuously cause errors
#   --jobs=200 # run up to 200 jobs at a time, does not limit number of cores used at a time
#   --cores=2000 # limit the number of cores used at a time to 2000, the limit for the 'short' partition on the cluster. Experimental.

#--- Finish up -----------------------------------------------------------------

set +x

end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo >&2 "Job ended at $end_time"

# Calculate and print the total runtime in hours and minutes
start_seconds=$(date -d "$start_time" +"%s")
end_seconds=$(date -d "$end_time" +"%s")
runtime_seconds=$((end_seconds - start_seconds))

runtime_hours=$((runtime_seconds / 3600))
runtime_minutes=$(((runtime_seconds % 3600) / 60))

echo >&2 "Total runtime: $runtime_hours hours and $runtime_minutes minutes"
echo >&2 "All done. You may have leave the cnakepit now. Have a nice day."
