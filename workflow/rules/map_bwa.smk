import os

rule get_ref:
    input:
        config['reference'+'_'+config['genome_version']] ["fasta_link"]
    output:
        config['reference'+'_'+config['genome_version']]["fasta"]
    shell:
        "mv {input} {output}"

rule get_ref_index:
    input:
        config['reference'+'_'+config['genome_version']] ["index_link"]
    output:
        config['reference'+'_'+config['genome_version']]["index"]
    shell:
        "mv {input} {output}"


if config["aligner"] == "bwa":

    rule bwa_index_reference:
        input:
            ref = ref_file,
        output:
            idx = multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        benchmark: "benchmarks/bwa_index.txt"
        log:
            "logs/bwa_index/bwa_index.log",
        params:
            algorithm="bwtsw",
        resources:
            mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
            slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
            runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
            cores=lambda wc, threads: threads
        log:
        conda:
            env_prefix + "trim_map" + env_suffix
        wrapper:
            "v1.7.0/bio/bwa/index"

    rule bwa_mem_samples:
        input:
            ref_index=config['reference'+'_'+config['genome_version']]["index"],
            reads=["results/trimmed/{sample}_1P.fq.gz", "results/trimmed/{sample}_2P.fq.gz"],
            idx=multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            temp("results/mapped/{sample}.bam") if not config["amplicon"] else "results/mapped/{sample}.bam"
        benchmark: 'benchmarks/bwa_mem/{sample}.txt'
        log:
            "logs/bwa_mem/{sample}.log",
        params:
            extra=r"-a -R '@RG\tID:{sample}\tSM:{sample}'", # -a to mark secondary alignments (for CNVkit)
            sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
            sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
            sort_extra="-@ {snakemake.threads}",  # Extra args for samtools/picard.
        threads: 16
        resources:
            mem=lambda wildcards, attempt: '%dG' % (12 * attempt),
            slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
            runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
            cores=lambda wc, threads: threads
        conda:
            env_prefix + "trim_map" + env_suffix
        wrapper:
            "v1.7.0/bio/bwa/mem"

elif config["aligner"] == "bwa-mem2":

    rule bwa_mem2_index_reference:
        input:
            ref = ref_file,
        output:
            idx = multiext(stem, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        benchmark: "benchmarks/bwa_mem2_index.txt"
        resources:
            mem=lambda wildcards, attempt: '%dG' % (80 * attempt),
            slurm_partition = lambda wildcards, attempt: 'highmem' if attempt > 1 else 'medium',
            runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
            cores=lambda wc, threads: threads
        log:
            "logs/bwa_mem2_index/bwa_index.log",
        conda:
            env_prefix + "trim_map" + env_suffix
        wrapper:
            "v5.0.2/bio/bwa-mem2/index"

    rule bwa_mem2_samples:
        input:
            ref_index=config['reference'+'_'+config['genome_version']]["index"],
            reads=["results/trimmed/{sample}_1P.fq.gz", "results/trimmed/{sample}_2P.fq.gz"],
            idx = multiext(stem, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        output:
            temp("results/mapped/{sample}.bam") if not config["amplicon"] else "results/mapped/{sample}.bam"
        benchmark: 'benchmarks/bwa_mem2/{sample}.txt'
        log: "logs/bwa_mem2/{sample}.log",
        params:
            extra=r"-a -R '@RG\tID:{sample}\tSM:{sample}'", # -a to mark secondary alignments (for CNVkit)
            stem_idx = stem,
        threads: 16
        resources:
            mem=lambda wildcards, attempt: '%dG' % (40 * attempt),
            slurm_partition = lambda wildcards, attempt: 'short' if attempt < 1 else 'medium',
            runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
            cores=lambda wc, threads: threads
        conda:
            env_prefix + "trim_map" + env_suffix
        shell:
            "bwa-mem2 mem "
            " -t {threads} "
            " {params.extra} "
            " {params.stem_idx} "
            " {input.reads} "
            " | samtools sort -@ {threads} -o {output} 2> {log}"


rule index_bwa_samples:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.bam.bai"
    benchmark: "benchmarks/samtools_index_bwa/{sample}.txt"
    log:
        "logs/samtools/index_bwa/{sample}.log"
    threads: 4
    resources:
        mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
    conda:
        env_prefix + "trim_map" + env_suffix
    params:
        extra="",  # optional params string
    wrapper:
        "v2.0.0/bio/samtools/index"

# if the data is based on hybrid-capture, flag PCR duplicates
if not config["amplicon"]:
    rule sambamba_mark_duplicates:
        input:
            BAMs_no_PCR_flags
        output:
            "results/mapped_marked/{sample}.bam",
            "results/mapped_marked/{sample}.bam.bai"
        benchmark: "benchmarks/sambamba_mark_duplicates/{sample}.txt"
        log:
            "logs/sambamba_mark_duplicates/{sample}.log"
        threads: 8
        resources:
            mem=lambda wildcards, attempt: '%dG' % (4 * attempt),
            runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
            slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
        conda:
            env_prefix + "cnv_calling" + env_suffix
        params:
            extra=""
        wrapper:
            "v2.0.0/bio/sambamba/markdup"

# merge multiple lanes of the same sample
rule merge_bam:
    input:
        lane1=lambda wildcards: 'results/mapped_marked/' + new_samples.at[wildcards.sample_stripped, 'lane1'] + '.bam',
        lane2=lambda wildcards: 'results/mapped_marked/' + new_samples.at[wildcards.sample_stripped, 'lane2'] + '.bam',
    output: "results/merged_bam/{sample_stripped}.bam"
    benchmark: "benchmarks/merge_bam/{sample_stripped}.txt"
    log: "logs/merge_bam/{sample_stripped}.log"
    resources: mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
    threads: 4
    conda: env_prefix + "trim_map" + env_suffix
    shell:
        "samtools merge --threads {threads} -o {output} {input.lane1} {input.lane2} 2> {log}"

# index the merged BAM files
rule index_merged_bam:
    input: "results/merged_bam/{sample_stripped}.bam"
    output: "results/merged_bam/{sample_stripped}.bam.bai"
    benchmark: "benchmarks/samtools_index_merged_bam/{sample_stripped}.txt"
    log: "logs/samtools/index_merged_bam/{sample_stripped}.log"
    threads: 4
    resources: mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
    conda: env_prefix + "trim_map" + env_suffix
    params: extra="",  # optional params string
    wrapper: "v2.0.0/bio/samtools/index"