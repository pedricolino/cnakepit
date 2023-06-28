# qualimap needs sorted bam files as input
rule sort_bwa:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/bam_sorted_bwa/{sample}_sorted.bam"
    log:
        "logs/samtools/sort_bwa/{sample}.log"
    threads:
        8
    conda:
        "../envs/qc_map.yaml"
    shell:
        "samtools sort -o {output} {input} -@ {threads} 2> {log}"

rule sort_bowtie2:
    input:
        "results/bowtie_mapped/{sample}.bam" 
    output:
        "results/bam_sorted_bowtie2/{sample}_sorted.bam"
    log:
        "logs/samtools/sort_bowtie2/{sample}.log"
    threads:
        8
    conda:
        "../envs/qc_map.yaml"
    shell:
        "samtools sort -o {output} {input} -@ {threads} 2> {log}"

rule qualimap_bwa:
    input:
       "results/bam_sorted_bwa/{sample}_sorted.bam"
    output:
        "results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html"
    log:
        "logs/qualimap_bwa/{sample}.log"
    params:
        outdir = "results/qc/qualimap/qualimap_bwa/{sample}"
    conda:
        "../envs/qc_map.yaml"
    shell:
        "qualimap bamqc -bam {input} -outdir {params.outdir} 2> {log}"

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
        "../envs/qc_map.yaml"
    shell:
        "qualimap bamqc -bam {input} -outdir {params.outdir}"


rule multiqc_map_bwa:
    input:
        expand("results/qc/qualimap/qualimap_bwa/{sample}/qualimapReport.html", sample = list(samples.index)),
    output:
        "results/qc_map_bwa/multiqc_report.html"
    log:
        "logs/qc_map_bwa/multiqc.log"
#    conda:
#        "../envs/qc_map.yaml"
#    shell:
#        "multiqc ./results/qc/qualimap/qualimap_bwa -o results/qc_map_bwa 2> {log}"
    wrapper:
      "v2.0.0/bio/multiqc"

rule multiqc_map_bowtie2:
    input:
        bowtie2 = expand("results/qc/qualimap/qualimap_bowtie2/{sample}/qualimapReport.html", sample = list(samples.index))
    output:
        "results/qc_map_bowtie2/multiqc_report.html"
    log:
        "logs/qc_map_bowtie2/multiqc.log"
    conda:
        "../envs/qc_map.yaml"
    shell:
        "multiqc ./results/qc/qualimap/qualimap_bowtie2 -o results/qc_map_bowtie2 2> {log}"

