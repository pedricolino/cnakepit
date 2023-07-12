rule purecn:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/{sample}.cnr',
        seg='results/cnvkit/{sample}.seg'
    output:
        "results/purecn/{sample}.rds"
    threads: 16
    log:
        "logs/purecn/{sample}.log",
    conda:
        "../envs/purecn.smk"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation Hclust"
    shell:
        "export PURECN='/fast/work/users/cemo10_c/miniconda/envs/purecn/lib/R/library/PureCN/extdata &&"
        "PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out {output} {params.extra}"