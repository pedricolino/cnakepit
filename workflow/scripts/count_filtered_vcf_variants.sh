#!/bin/bash

output_file="results/stats/mutect2_vcf_filtered_variants_counts.tsv"

# Create the TSV file and add the header
echo -e "SampleName\tTotalVariants\tPassedVariants\tGermlineVariants" > "$output_file"

# Get the total number of files
total_files=$(ls results/mutect2/filtered/*.vcf.gz | wc -l)
processed_files=0

for file in results/mutect2/filtered/*.vcf.gz; do
    # Count the total number of variants
    total_variants=$(zcat "$file" | grep '^[^#]' | wc -l)

    # Count the number of variants with PASS in the filter column
    pass_variants=$(zcat "$file" | grep '^[^#]*PASS' | wc -l)

    # Count the number of variants with germline in the filter column
    germline_variants=$(zcat "$file" | grep '^[^#]*germline' | wc -l)

    SampleName=$(basename "$file" | sed -e 's/_filtered.vcf.gz//g')

    # Append the counts to the TSV file
    echo -e "$SampleName\t$total_variants\t$pass_variants\t$germline_variants" >> "$output_file"

    # Update the progress
    processed_files=$((processed_files + 1))
    echo -ne "Progress: $processed_files/$total_files\r"
done

echo -e "\nProcessing completed!"
