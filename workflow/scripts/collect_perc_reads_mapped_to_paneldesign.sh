#!/bin/bash

output_file="results/stats/mapped_to_reference.txt"
input_folder="results/qc/qualimap/qualimap_bwa"

touch $output_file

for file in $(ls -R $input_folder/*/*.html)
do
    cat $file \
        | sed -n '181p' \
        | sed -e "s/.*  \///g" \
        | sed -e "s/%<\/td>//g" \
        >> $output_file
done

echo "File written to $output_file"