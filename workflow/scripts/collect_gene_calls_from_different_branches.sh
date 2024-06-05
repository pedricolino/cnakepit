#!/bin/bash

# Search for all branches in results folder. Every folder starts with purecn, all branches are its subfolders.
branches=$(ls -d results/*purecn*/*/)

# default output folder
OUT=results/stats/gene_calls

# check whether the output folder is writable, otherwise change it
if [ -w results/stats ]; then 
    echo "Default output path $OUT is writable and will be used."
else 
    echo "Default output path $OUT is not writable, changing to results_all_write_permissions/stats/gene_calls"
    OUT=results_all_write_permissions/stats/gene_calls
fi

# Create the output folder if it does not exist
mkdir -p $OUT

# Iterate over each branch and call collect_gene_calls.sh
for branch in $branches; do
    # Extract the branch name by removing the results/ prefix and the trailing slash, and replacing / with _
    branch_name=$(echo $branch | sed 's/results\///g' | sed 's/\//_/g' | sed 's/_$//g')
    # Call collect_gene_calls.sh with the branch name
    bash collect_gene_calls.sh -o ${OUT}/all_gene_calls_${branch_name}_gene_calls.csv -d $branch -s resources/data/panel_samples_with_wes_counterparts.csv
done
