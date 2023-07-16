rule purecn:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/{sample}.cnr',
        seg='results/cnvkit/{sample}.seg',
        #script="workflow/scripts/purecn.R"
    output:
        "results/purecn/{sample}/{sample}_loh.csv"
        "results/purecn/{sample}/{sample}_chromosomes.pdf"
        "results/purecn/{sample}/{sample}_variants.csv"
        "results/purecn/{sample}/{sample}_local_optima.pdf"
        "results/purecn/{sample}/{sample}_genes.csv"
        "results/purecn/{sample}/{sample}_dnacopy.seg"
        "results/purecn/{sample}/{sample}.pdf"
        "results/purecn/{sample}/{sample}_segmentation.pdf"
        "results/purecn/{sample}/{sample}.rds"
        "results/purecn/{sample}/{sample}.csv"
    threads: 16
    log:
        "logs/purecn/{sample}.log",
    conda:
        "../envs/purecn.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation Hclust"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out {output} {params.extra}
        """