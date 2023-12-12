#!/bin/bash

dir=$1
out=$2

# Print the header from the first .icns file
echo -e "chromosome\tstart\tend\tgene\tlog2\tbaf\tcn\tcn1\tcn2\tdepth\tprobes\tweight\tsample" > $out

# Process each .icns file
total_files=$(ls -1 "$dir"/*.icns | wc -l)
current_file=1

for file in "$dir"/*.icns; do
    # extract sample name
    name=$(basename "$file" | sed -e 's/\.icns//g')

    # Skip the header and append the filename to each row
    tail -n +2 "$file" | awk -v OFS="\t" -v file="$name" '{print $0, file}' >> $out

    # Print progress
    echo -ne "Processing file $current_file of $total_files\r"
    ((current_file++))
done

echo -e "\nProcessing completed!"