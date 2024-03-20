#--- Downloads start ---#

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule download_gnomad:
    input: HTTP.remote(config['gnomad_af_only'+'_'+config['genome_version']]["vcf_link"], keep_local=True)
    output: config['gnomad_af_only'+'_'+config['genome_version']]["vcf"]
    shell: "mv {input} {output}"

rule download_gnomad_index:
    input: HTTP.remote(config['gnomad_af_only'+'_'+config['genome_version']]["index_link"], keep_local=True)
    output: config['gnomad_af_only'+'_'+config['genome_version']]["index"]
    shell: "mv {input} {output}"

rule download_common_biallelic:
    input: HTTP.remote(config['common_germline_variants'+'_'+config['genome_version']]["vcf_link"], keep_local=True)
    output: config['common_germline_variants'+'_'+config['genome_version']]["vcf"]
    shell: "mv {input} {output}"

rule download_common_biallelic_index:
    input: HTTP.remote(config['common_germline_variants'+'_'+config['genome_version']]["index_link"], keep_local=True)
    output: config['common_germline_variants'+'_'+config['genome_version']]["index"]
    shell: "mv {input} {output}"

if config['compute_dict']:
    rule mutect2_ref_dict:
        input: config['reference'+'_'+config['genome_version']]["fasta"]
        output: config['reference'+'_'+config['genome_version']]["dict"]
        conda: "../envs/cnv_calling.yaml"
        log: "logs/mutect2_dict/mutect2_dict.log"
        shell: "gatk CreateSequenceDictionary -R {input} 2> {log}"
else:
    # Get the dict from the web. If ref. is masked, the .dict changes but it is not used anymore after mapping.
    rule download_ref_dict:
        input: HTTP.remote(config['reference'+'_'+config['genome_version']]["dict_link"], keep_local=True)
        output: config['reference'+'_'+config['genome_version']]["dict"]
        benchmark: "benchmarks/download_ref_dict.log"
        shell: "mv {input} {output}"

#--- Downloads end ---#

rule mutect2_bam:
    input:
        fasta=config['reference'+'_'+config['genome_version']]["fasta"],
        map=BAMs_for_CNV_calling,
        idx=BAM_index_for_CNV_calling,
        dict=config['reference'+'_'+config['genome_version']]["dict"],
        targets=config["panel_design"],
        gnomad=config['gnomad_af_only'+'_'+config['genome_version']]["vcf"],
        gnomad_index=config['gnomad_af_only'+'_'+config['genome_version']]["index"]
    output:
        vcf="results/mutect2/unfiltered/{sample}.vcf.gz",
        f1r2="results/mutect2/f1r2/{sample}.tar.gz"
    benchmark: "benchmarks/mutect2_bam/{sample}.txt"
    priority: -1 # Mutect2 is very slow, so we want to run downstream rules of already Mutect2-processed samples first
    params:
        extra="--genotype-germline-sites true --genotype-pon-sites true --interval-padding 75"
    threads: 16 # Confirmed in log files that it works, see https://www.biostars.org/p/9549710/#9550707
    resources:
        mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
        slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
        runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
        cores=lambda wc, threads: threads
    log:
        "logs/mutect2/{sample}.log",
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.map} -L {input.targets} -O {output.vcf} {params.extra} --native-pair-hmm-threads {threads} --germline-resource {input.gnomad} --f1r2-tar-gz {output.f1r2} --tmp-dir ${{TMPDIR}} &> {log}"

rule learn_read_orientation_model:
    input:
        f1r2="results/mutect2/f1r2/{sample}.tar.gz",
        idx=BAM_index_for_CNV_calling
    output:
        "results/mutect2/read_orientation_model/{sample}.tar.gz"
    resources:
        mem=lambda wildcards, attempt: '%dG' % (2 * attempt),
    benchmark: "benchmarks/read_orientation_model/{sample}.txt"
    log:
        "logs/read_orientation_model/{sample}.log"
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        "gatk LearnReadOrientationModel -I {input.f1r2} -O {output} --tmp-dir ${{TMPDIR}} &> {log}"

rule get_pile_up_summaries:
    input:
        bam=BAMs_for_CNV_calling,
        bam_idx=BAM_index_for_CNV_calling,
        common=config['common_germline_variants'+'_'+config['genome_version']]["vcf"],
        common_index=config['common_germline_variants'+'_'+config['genome_version']]["index"],
    output:
        "results/mutect2/pile_up_summaries/{sample}.table"
    benchmark: "benchmarks/pile_up_summaries/{sample}.txt"
    log:
        "logs/pile_up_summaries/{sample}.log"
    conda:  
        "../envs/cnv_calling.yaml"
    shell:
        "gatk GetPileupSummaries -I {input.bam} -V {input.common} -L {input.common} -O {output} --tmp-dir ${{TMPDIR}} &> {log}"

rule calculate_contamination:
    input:
        pileup="results/mutect2/pile_up_summaries/{sample}.table",
        reference=config['reference'+'_'+config['genome_version']]["fasta"],
        dict=config['reference'+'_'+config['genome_version']]["dict"],
    output:
        "results/mutect2/contamination/{sample}.txt"
    benchmark: "benchmarks/calculate_contamination/{sample}.txt"
    log:
        "logs/calculate_contamination/{sample}.log"
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        "gatk CalculateContamination -I {input.pileup} -O {output} --tmp-dir ${{TMPDIR}} &> {log}"

rule filter_mutect_calls:
    input:
        vcf="results/mutect2/unfiltered/{sample}.vcf.gz",
        reference=config['reference'+'_'+config['genome_version']]["fasta"],
        rom="results/mutect2/read_orientation_model/{sample}.tar.gz",
        contamination_table="results/mutect2/contamination/{sample}.txt"
    log:
        "logs/filter_mutect_calls/{sample}.log",
    conda:
        "../envs/cnv_calling.yaml"
    output:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
    benchmark: "benchmarks/filter_mutect_calls/{sample}.txt"
    shell:
        "gatk FilterMutectCalls -R {input.reference} -V {input.vcf} --orientation-bias-artifact-priors {input.rom} --contamination-table {input.contamination_table}  -O {output.vcf_filt} --tmp-dir ${{TMPDIR}} &> {log}"

# for CNVkit, extract only variants classified as germline
rule extract_germline_variants:
    input:
        "results/mutect2/filtered/{sample}_filtered.vcf.gz"
    log: 
        "logs/extract_germline_variants/{sample}.log"
    output:
        "results/mutect2/germline/{sample}_germline.vcf"
    benchmark:
        "benchmarks/extract_germline_variants/{sample}.txt"
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        """
        bcftools view -i 'FILTER~"germline"' {input} > {output} 2> {log}
        """