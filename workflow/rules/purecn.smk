# required for patched PSCBS segm. method
rule install_lima1_pscbs:
    output:
        "results/purecn/pscbs_check"
    conda:
        "../envs/primary_env.yaml"
    log:
        "logs/purecn/install_lima1_pscbs.log"
    shell:
        "Rscript workflow/scripts/install_lima1_pscbs.R && touch {output} 2> {log}"


########## Segmentation with different methods ##########

rule purecn_cbs_pscbs:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/general/{sample}.cnr',
        seg='results/cnvkit/cbs/{sample}.seg',
        install="results/purecn/pscbs_check"
    output:
        "results/purecn/cbs_pscbs/{sample}/{sample}.rds"
    benchmark: "benchmarks/purecn/cbs_pscbs/{sample}.txt"
    log:
        "logs/purecn/cbs_pscbs/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation PSCBS --sex F"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out results/purecn/cbs_pscbs/{params.sampleid}/{params.sampleid} {params.extra} --genome hg38 &> {log}
        """

# "The --stats-file is only supported for Mutect 1.1.7. Mutect 2 provides the filter flags directly in the VCF.""
# The --fun-segmentation argument controls if the data should to be re-segmented using germline BAFs (default). Set this value to none if the provided segmentation should be used as is. The recommended Hclust will only cluster provided segments.

rule purecn_hmm_pscbs:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/general/{sample}.cnr',
        seg='results/cnvkit/hmm/{sample}.seg',
        install="results/purecn/pscbs_check"
    output:
        "results/purecn/hmm_pscbs/{sample}/{sample}.rds"
    log:
        "logs/purecn/hmm_pscbs/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation PSCBS --sex F"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out results/purecn/hmm_pscbs/{params.sampleid}/{params.sampleid} {params.extra} --genome hg38 &> {log}
        """

rule purecn_cbs_hclust:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/general/{sample}.cnr',
        seg='results/cnvkit/cbs/{sample}.seg',
        #script="workflow/scripts/purecn.R"
    output:
        "results/purecn/cbs_hclust/{sample}/{sample}.rds"
    benchmark: "benchmarks/purecn/cbs_hclust/{sample}.txt"
    log:
        "logs/purecn/cbs_hclust/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation Hclust --sex F"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out results/purecn/cbs_hclust/{params.sampleid}/{params.sampleid} {params.extra} --genome hg38 &> {log}
        """

rule purecn_hmm_hclust:
    input:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
        copy_ratios='results/cnvkit/general/{sample}.cnr',
        seg='results/cnvkit/hmm/{sample}.seg',
    output:
        "results/purecn/hmm_hclust/{sample}/{sample}.rds"
    log:
        "logs/purecn/hmm_hclust/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    params:
        sampleid="{sample}",
        extra="--genome hg38 --force --postoptimize --seed 123 --funsegmentation Hclust --sex F"
    shell:
        """PURECN=$(Rscript -e "cat(system.file('extdata', package = 'PureCN'))")
        Rscript $PURECN/PureCN.R --vcf {input.vcf_filt} --sampleid {params.sampleid} --tumor {input.copy_ratios} --segfile {input.seg} --out results/purecn/hmm_hclust/{params.sampleid}/{params.sampleid} {params.extra} --genome hg38 &> {log}
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
