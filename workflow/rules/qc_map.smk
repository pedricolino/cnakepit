################ BWA ######################

# qualimap needs sorted bam files as input
rule sort_bwa:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/bam_sorted_bwa/{sample}_sorted.bam"
    benchmark:
        "benchmarks/sort_bwa/{sample}.txt"
    log:
        "logs/samtools/sort_bwa/{sample}.log"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
        runtime=24*60, # 24h
        slurm_partition='medium'
    conda:
        "../envs/primary_env.yaml"
    shell:
        "samtools sort -o {output} {input} -@ {threads} 2> {log}" 

rule qualimap_bwa:
    input:
       "results/bam_sorted_bwa/{sample}_sorted.bam"
    output:
        "results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html"
    benchmark:
        "benchmarks/qualimap_bwa/{sample}.txt"
    threads: 16
    priority: -2 # error-prone rule, run it last
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=5*24*60, # 24h seems to be not enough for many samples
        slurm_partition='medium'
    log:
        "logs/qualimap_bwa/{sample}.log"
    params:
        outdir = "results/qc/qualimap/qualimap_bwa/{sample}"
    conda:
        "../envs/qc_map.yaml"
    shell:
        "qualimap bamqc -bam {input} -nt {threads} --java-mem-size={resources.mem} -outdir {params.outdir} 2> {log}"

rule multiqc_map_bwa:
    input:
        expand("results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html", sample = list(samples.index))
    output:
        "results/qc_map_bwa/multiqc_report.html"
    benchmark:
        "benchmarks/multiqc_map_bwa.txt"
    priority: -2 # error-prone rule, run it last
    resources:
        mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
        runtime=24*60, # 24h
        slurm_partition='medium'
    log:
        "logs/qc_map_bwa/multiqc.log"
    conda:
        "../envs/primary_env.yaml"
    shell:
        "multiqc ./results/qc/qualimap/qualimap_bwa -o results/qc_map_bwa 2> {log}"

################# Bowtie2 ######################

rule sort_bowtie2:
    input:
        "results/bowtie_mapped/{sample}.bam" 
    output:
        "results/bam_sorted_bowtie2/{sample}_sorted.bam"
    log:
        "logs/samtools/sort_bowtie2/{sample}.log"
    threads:
        16
    conda:
        "../envs/primary_env.yaml"
    shell:
        "samtools sort -o {output} {input} -@ {threads} 2> {log}"  

rule qualimap_bowtie2:
    input:
        "results/bam_sorted_bowtie2/{sample}_sorted.bam"
    output:
        "results/qc/qualimap/qualimap_bowtie2/{sample}/qualimapReport.html"
    log:
        "logs/qualimap_bowtie2/{sample}.log"
    params:
        outdir = "results/qc/qualimap/qualimap_bowtie2/{sample}"
    conda:
        "../envs/primary_env.yaml"
    shell:
        "qualimap bamqc -bam {input} -outdir {params.outdir}"

rule multiqc_map_bowtie2:
    input:
        bowtie2 = expand("results/qc/qualimap/qualimap_bowtie2/{sample}/qualimapReport.html", sample = list(samples.index))
    output:
        "results/qc_map_bowtie2/multiqc_report.html"
    log:
        "logs/qc_map_bowtie2/multiqc.log"
    conda:
        "../envs/primary_env.yaml"
    shell:
        "multiqc ./results/qc/qualimap/qualimap_bowtie2 -o results/qc_map_bowtie2 2> {log}"

