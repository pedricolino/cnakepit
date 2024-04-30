#!/bin/bash

sample_names_table=$1
OUT=results/collection_of_specific_samples
mkdir -p $OUT

# for every folder starting with results/purecn
for pon in results/purecn*; do
    echo "Processing $pon"
    pon_name=$(basename $pon)
    mkdir -p $OUT/$pon_name
    for branch in $(ls -d $pon/*/); do
        branch_name=$(basename $branch)
        mkdir -p $OUT/$pon_name/$branch_name
        while IFS= read -r name; do
            seg_files=$(find $branch -name "$name*_dnacopy.seg")
            if [ -n "$seg_files" ]; then
                cp "$seg_files" "$OUT/$pon_name/$branch_name/${name}_dnacopy.seg"
            fi

            genes_files=$(find $branch -name "$name*_genes.csv")
            if [ -n "$genes_files" ]; then
                cp "$genes_files" "$OUT/$pon_name/$branch_name/${name}_genes.csv"
            fi

            rds_files=$(find $branch -name "$name*.rds")
            if [ -n "$rds_files" ]; then
                cp "$rds_files" "$OUT/$pon_name/$branch_name/${name}.rds"
            fi
        done <$sample_names_table
    done
done
