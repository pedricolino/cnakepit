rule index:
    input:
        "results/bam_sorted_bwa/{sample}_sorted.bam"
    output:
        "results/bam_sorted_bwa/{sample}_sorted.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    threads:
        4
    conda:
        "../envs/stats.yaml"
    shell:
        "samtools index {input} -@ {threads} 2> {log}"

rule stats:
    input:
        "results/bam_sorted_bwa/{sample}_sorted.bam",
        "results/bam_sorted_bwa/{sample}_sorted.bam.bai"
    output:
        "results/stats/{sample}.stats"
    log:
        "logs/samtools/stats/{sample}.log"
    conda:
        "../envs/stats.yaml"
    shell:
        "samtools idxstats {input[0]} > {output} 2> {log}"


rule rpk:
    input:
        "results/stats/{sample}.stats"
    output:
        "results/stats/{sample}.stats_aug"
    script:
        "../scripts/calc_rpk.py"



