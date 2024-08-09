rule stats:
    input:
        "results/bam_sorted_bwa/{sample}_sorted.bam",
        "results/bam_sorted_bwa/{sample}_sorted.bam.bai"
    output:
        "results/stats/{sample}.stats"
    benchmark: "benchmarks/stats/{sample}.stats.benchmark.txt"
    log:
        "logs/samtools/stats/{sample}.log"
    conda:
        env_prefix + 'qc' + env_suffix
    shell:
        "samtools idxstats {input[0]} > {output} 2> {log}"


rule rpk:
    input:
        "results/stats/{sample}.stats"
    output:
        "results/stats/{sample}.stats_aug"
    benchmark: "benchmarks/stats/{sample}.stats_aug.benchmark.txt"
    script:
        "../scripts/calc_rpk.py"



