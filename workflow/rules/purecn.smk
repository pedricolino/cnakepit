rule purecn:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/{sample}.cnr',
        seg='results/cnvkit/{sample}.seg',
        #script="workflow/scripts/purecn.R"
    output:
        "results/purecn/{sample}/{sample}.rds"
    threads: 16
    log:
        "logs/purecn/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation Hclust"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out results/purecn/{params.sampleid} {params.extra} && mv results/purecn/{params.sampleid}/{params.sampleid}.log logs/purecn/{params.sampleid}.log
        """

# rule check_end:
#     input:
#         expand('results/purecn/{sample}/{sample}.rds', sample=samples.index)
#     output:
#         "results/.check_end"
#     log:
#         "logs/purecn/check_end.log"
#     shell:
#         "touch {output} 2> {log}"