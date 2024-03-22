#!/bin/bash

mkdir -p results/stats

printf "\nCollect subset of specific (e.g. WES-matching) samples in one folder (results/collection_of_specific_samples).\n"
./workflow/scripts/collect_purecn_results_for_specific_samples.sh ~/work/sign-oc/wes_panel_results/corresponding_panels_names.tsv

printf "\nCount the number of successful PureCN runs\n"
./workflow/scripts/count_successful_purecn_runs.sh

printf "\nCollect for all samples the percentage of reads mapped to the panel design\n"
./workflow/scripts/collect_perc_reads_mapped_to_paneldesign.sh

printf "\nCompute mean gene coverage for all samples with CNVkit's CNN output files\n"
Rscript workflow/scripts/calc_mean_gene_coverage.R

printf "\nCollect all benchmarks into one file\n"
python workflow/scripts/collect_benchmarks.py

printf "\nCount the total number of variants and the number of variants with PASS in the filter column\n"
./workflow/scripts/count_filtered_vcf_variants.sh

printf "\nCount SNVs found by Mutect2\n"
./workflow/scripts/count_snv.sh

printf "\nCollect all PureCN optima for tumor purity and ploidy information for all samples\n"
Rscript workflow/scripts/collect_purecn_optima.R $(ls -d results/purecn/*/ | xargs)

printf "\nCollect for all genes and samples the absolute copy numbers and log2 ratios found by PureCN\n"
Rscript workflow/scripts/collect_gene_copy_numbers_log2ratios.R $(ls -d results/purecn/*/ | xargs)

