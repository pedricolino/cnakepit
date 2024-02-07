#!/bin/bash

# Specify the directory
dir="results/mutect2/filtered"

# Specify the output file
output_file="results/stats/mutect2_snv_counts.tsv"

echo -e "sample\tref\talt\tcount" > "$output_file"

# Get the total number of files in the directory
total_files=$(ls -1 "$dir"/*.vcf.gz | wc -l)

# Initialize a counter for processed files
processed_files=0

# Loop over all VCF files in the directory
for file in "$dir"/*.vcf.gz; do
    # Increment the counter for processed files
    ((processed_files++))

    # Print the progress
    echo -ne "Processing file $processed_files/$total_files\r"

    # Loop over all possible SNVs
    for ref in C G T A; do
        for alt in C G T A; do
            # Skip if ref and alt are the same
            if [ "$ref" != "$alt" ]; then
                # Count the number of occurrences of the SNV
                count=$(zcat "$file" | grep -v '^#' | awk -v ref="$ref" -v alt="$alt" '$4 == ref && $5 == alt' | wc -l)
                
                # extract sample name
                name=$(basename "$file" | sed -e 's/_filtered.vcf.gz//g')

                # Write the SNV and the count to the output file
                echo -e "$name\t$ref\t$alt\t$count" >> "$output_file"
            fi
        done
    done
done

echo -e "\nProcessing completed!"
