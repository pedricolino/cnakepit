import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_ref:
    input:
        HTTP.remote(config['reference'+'_'+config['genome_version']] ["fasta_link"], keep_local=True)
    output:
        config['reference'+'_'+config['genome_version']]["fasta"]
    shell:
        "mv {input} {output}"

rule get_ref_index:
    input:
        HTTP.remote(config['reference'+'_'+config['genome_version']] ["index_link"], keep_local=True)
    output:
        config['reference'+'_'+config['genome_version']]["index"]
    shell:
        "mv {input} {output}"

rule bwa_index_reference:
    input:
        ref = ref_file,
    output:
        idx=multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    benchmark: "benchmarks/bwa_index.txt"
    log:
        "logs/bwa_index/bwa_index.log",
    params:
        algorithm="bwtsw",
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

rule bwa_index_samples:
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
            mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
        conda:
            env_prefix + "cnv_calling" + env_suffix
        params:
            extra=""
        wrapper:
            "v2.0.0/bio/sambamba/markdup"
