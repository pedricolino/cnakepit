#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate R

printf "\n1. Collect subset of specific (e.g. WES-matching) samples in one folder (results/collection_of_specific_samples).\n"
./workflow/scripts/collect_purecn_results_for_specific_samples.sh resources/data/names_of_subset_samples.csv

mkdir -p results/stats

printf "\n2. Collect the gene-level copy number changes for all samples per branch.\n"
Rscript workflow/scripts/gene_calls_per_branch_to_tsv.R

printf "\n3. Count the number of successful PureCN runs\n"
./workflow/scripts/count_successful_purecn_runs.sh

printf "\n4. Collect for all samples the percentage of reads mapped to the panel design\n"
./workflow/scripts/collect_perc_reads_mapped_to_paneldesign.sh

printf "\n5. Compute mean gene coverage for all samples with CNVkit's CNN output files\n"
Rscript workflow/scripts/calc_mean_gene_coverage.R

printf "\n6. Collect all benchmarks into one file\n"
Rscript workflow/scripts/collect_benchmarks.R

printf "\n7. Render the benchmarks report\n"
Rscript -e "rmarkdown::render('workflow/scripts/summarise_benchmarks.Rmd',  output_file = 'results/stats/benchmarks_summary.html')"

printf "\n8. Count the total number of variants and the number of variants with PASS in the filter column\n"
./workflow/scripts/count_filtered_vcf_variants.sh

printf "\n9. Count SNVs found by Mutect2\n"
./workflow/scripts/count_snv.sh

printf "\n10. Collect all PureCN optima for tumor purity and ploidy information for all samples\n"
Rscript workflow/scripts/collect_purecn_optima.R $(ls -d results/purecn/*/ | xargs)

printf "\n11. Collect for all genes and samples the absolute copy numbers and log2 ratios found by PureCN\n"
Rscript workflow/scripts/collect_gene_copy_numbers_log2ratios.R $(ls -d results/purecn/*/ | xargs)

cp results/qc_map_bwa/multiqc_data/multiqc_general_stats.txt results/stats/multiqc_general_stats.txt
cp results/qc_map_bwa/multiqc_data/general_stats_table.tsv results/stats/general_stats_table.tsv

ls results/purecn/cbs_Hclust >results/stats/all_processed_samples.csv
