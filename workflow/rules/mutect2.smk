#ref_ref = str(Path("results") / "reference" / ref_file.stem)

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_gnomad:
    input:
        HTTP.remote("www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz", keep_local=True)
    output:
        config["germline-resource"]
    shell:
        "mv {input} {output}"

rule get_gnomad_index:
    input:
        HTTP.remote("www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz.tbi", keep_local=True)
    output:
        config["germline-resource-index"]
    shell:
        "mv {input} {output}"

rule mutect2_bam:
    input:
        fasta=config["ref"],
        map="results/bam_sorted_bwa/{sample}_sorted.bam",
        dict="resources/reference/hg38.dict",
        idx="results/bam_sorted_bwa/{sample}_sorted.bam.bai",
        targets=config["bed_w_chr"],
        gnomad=config["germline-resource"],
        gnomad_index=config["germline-resource-index"]
    output:
        vcf="results/mutect2/variant_bam/{sample}.vcf.gz",
        #bam="results/mutect2/variant_bam/{sample}.bam",
    #message:
        #"Testing Mutect2 with {wildcards.sample}"
    params:
        extra="--genotype-germline-sites true --genotype-pon-sites true --interval-padding 75"
    threads: 16
    log:
        "logs/mutect2/{sample}.log",
    conda:
        "../envs/mutect2.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.map} -L {input.targets} -O {output.vcf} {params.extra} --germline-resource {input.gnomad} 2> {log}"

# rule mutect2_ref_dict:
#     input:
#         config["ref"]
#     output:
#         "resources/reference/hg38.dict"
#     conda:
#         "../envs/stats.yaml"
#     log:
#         "logs/mutect2_dict/mutect2_dict.log"
#     shell:
#         "samtools dict {input} > {output} 2> {log}"

rule get_ref_dict:
    input:
        HTTP.remote("storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict", keep_local=True)
    output:
        "resources/reference/hg38.dict"
    shell:
        "mv {input} {output}"

rule filter_mutect_calls:
    input:
        vcf="results/mutect2/variant_bam/{sample}.vcf.gz",
        reference=config["ref"],
    threads: 16
    log:
        "logs/filter_mutect_calls/{sample}.log",
    conda:
        "../envs/mutect2.yaml"
    output:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
    shell:
        "gatk FilterMutectCalls -R {input.reference} -V {input.vcf} -O {output.vcf_filt} 2> {log}"