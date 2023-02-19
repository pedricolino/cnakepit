run here with "snakemake use-conda --cores 4" to ensure the environment is being created
this is needed if run on a vanilla machine with snakemake only.
if no snakemake installed, create a new conda env first with "conda install -n <your just created basic env> -c conda-forge mamba"
"mamba create -c conda-forge -c bioconda -n <env name> snakemake"
"conda activate <env name>

locally: conda activate snakemake-vanilla
