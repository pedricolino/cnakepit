#!/bin/bash

# The medium project/queue is a sensible default.
#SBATCH --partition medium
# Set a required running time for the master job.
#SBATCH --time=7-00:00:00
# Reserve some resources
#SBATCH --mem=64G
#SBATCH --ntasks=96
# Keep current environment variables
#SBATCH --export=all
# Send a mail upon job completion and error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user cedric.moris@bih-charite.de
# Logs should be written into "slurm_log" sub directory.
#SBATCH --output logs/slurm_log/%x-%J.log
# Use more descriptive name in Slurm.
#SBATCH --job-name 5samples


# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Activate bash cmd printing, debug info ------------------------------------

set -x
>&2 hostname
>&2 date

# Activate conda environment ------------------------------------------------

# clutters log files and should test snakemake -n first anyway
#eval "$(conda shell.bash hook)"
#conda activate snakemake-vanilla

# Kick off Snakemake --------------------------------------------------------

snakemake --use-conda --cores 96

# Print date after finishing, for good measure ------------------------------

>&2 date
>&2 echo "All done. Have a nice day."
