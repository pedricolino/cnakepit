#!/bin/bash

sample_names_table=$1
OUT=results/collection_of_specific_samples
mkdir -p $OUT

# for every folder starting with results/purecn
for pon in results/purecn*; do
    echo "Processing $pon"
    for branch in $(ls -d $pon/*/); do
        mkdir -p $OUT/$branch
        while IFS= read -r name; do
            cp $(find $branch -name "$name*_dnacopy.seg") "results/collection_of_specific_samples/${branch}${name}_dnacopy.seg"
            cp $(find $branch -name "$name*_genes.csv") "results/collection_of_specific_samples/${branch}${name}_genes.csv"
            cp $(find $branch -name "$name*.rds") "results/collection_of_specific_samples/${branch}${name}.rds"
        done <$sample_names_table
    done
done
