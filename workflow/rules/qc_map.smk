################ BWA ######################

rule bed_for_qualimap:
    input: config["panel_design"]
    output: config["qualimap_bed"]
    shell:
        """
        cat {input} | grep '^chr' | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4,0,"."}}' > {output}
        """ # https://stackoverflow.com/questions/50965417/awk-command-fails-in-snakemake-use-singularity

# qualimap needs sorted bam files as input
rule qualimap_bwa:
    input:
       map=BAMs_for_CNV_calling,
        targets=config["qualimap_bed"]
    output:
        "results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html"
    benchmark:
        "benchmarks/qualimap_bwa/{sample}.txt"
    threads: 4
    priority: -2 # error-prone rule, run it last
    resources:
        mem=lambda wildcards, attempt: '%dG' % (12 * attempt),
        slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
        runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60,
        cores=lambda wc, threads: threads
    log:
        "logs/qualimap_bwa/{sample}.log"
    params:
        outdir = "results/qc/qualimap/qualimap_bwa/{sample}"
    conda:
        "../envs/qc.yaml"
    shell:
        "qualimap bamqc -bam {input.map} --feature-file {input.targets} -nt {threads} --java-mem-size={resources.mem} -outdir {params.outdir} 2> {log}"

rule multiqc_map_bwa:
    input:
        expand("results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html", sample = list(samples.index))
    output:
        "results/qc_map_bwa/multiqc_report.html"
    benchmark:
        "benchmarks/multiqc_map_bwa.txt"
    priority: -2 # error-prone rule, run it last
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        runtime=24*60, # 24h
        slurm_partition='medium'
    log:
        "logs/qc_map_bwa/multiqc.log"
    conda:
        "../envs/qc.yaml"
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
        "../envs/qc.yaml"
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
        "../envs/qc.yaml"
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
        "../envs/qc.yaml"
    shell:
        "multiqc ./results/qc/qualimap/qualimap_bowtie2 -o results/qc_map_bowtie2 2> {log}"

