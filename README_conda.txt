test
run here with "snakemake --use-conda --cores 8" to ensure the environment is being created
this is needed if run on a vanilla machine with snakemake only.

if no snakemake installed, create a new conda env first with "conda install -n <your just created basic env> -c conda-forge mamba"
"mamba create -c conda-forge -c bioconda -n <env name eg. snakemake-vanilla> snakemake"
"conda activate <env name>

on the cluster: only run while connected to computing node with "srun --time 4-00 --mem=8G --ntasks=8 --pty bash -i"

locally: conda activate snakemake-vanilla
