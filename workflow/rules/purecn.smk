from random import randint

# required for patched PSCBS segm. method
rule install_lima1_PSCBS:
    output: 'results/purecn/PSCBS_check'
    conda: env_prefix + 'cnv_calling' + env_suffix
    log: 'logs/purecn/install_lima1_PSCBS.log'
    shell: 'Rscript workflow/scripts/install_lima1_PSCBS.R && touch {output} 2> {log}'


# Segmentation with default methods CBS from CNVkit and Hclust from PureCN

rule purecn_cbs_Hclust:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'], # simpleRepeats track as in PureCN's vignette
        seg='results/cnvkit'+suffix+'/cbs/{sample}.seg'
    output: 'results/purecn'+suffix+'/cbs_Hclust/{sample}/{sample}.csv'
    resources: mem=lambda wildcards, attempt: '%dG' % (4 * 8 * attempt), # 4 GB per thread, 8 threads
    threads: 8
    benchmark: 'benchmarks/purecn'+suffix+'/cbs_Hclust/{sample}.txt'
    log: 'logs/purecn'+suffix+'/cbs_Hclust/{sample}.log',
    priority: 10 # run as early as possible to check results of first few samples
    conda: env_prefix + 'cnv_calling' + env_suffix
    params:
        cnvkit_method='cbs',
        purecn_method='Hclust',
        sampleid='{sample}',
        random_nb=randint(1,1000), # use a different seed on retry in case of failure
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']
    shell:
        '''PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
        Rscript $PURECN/PureCN.R \
            --vcf {input.vcf_filt} \
            --sampleid {params.sampleid} \
            --tumor {input.copy_ratios} \
            --seg-file {input.seg} \
            --snp-blacklist {input.blacklist} \
            --out results/purecn{params.suffix}/{params.cnvkit_method}_{params.purecn_method}/{params.sampleid}/{params.sampleid} \
            --genome {params.genome} \
            {params.sex} \
            --fun-segmentation {params.purecn_method} \
            --min-base-quality 20 \
            --seed {params.random_nb} \
            --force \
            --post-optimize \
            --cores {threads} \
            &> {log}
        '''

# 'The --stats-file is only supported for Mutect 1.1.7. Mutect 2 provides the filter flags directly in the VCF.''
# The --fun-segmentation argument controls if the data should to be re-segmented using germline BAFs (default). Set this value to none if the provided segmentation should be used as is. The recommended Hclust will only cluster provided segments.


# Segmentation with other methods

use rule purecn_cbs_Hclust as purecn_cbs_PSCBS with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
        seg='results/cnvkit'+suffix+'/cbs/{sample}.seg',
        install='results/purecn/PSCBS_check'
    output: 'results/purecn'+suffix+'/cbs_PSCBS/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/cbs_PSCBS/{sample}.txt'
    log: 'logs/purecn'+suffix+'/cbs_PSCBS/{sample}.log',
    params:
        cnvkit_method='cbs',
        purecn_method='PSCBS',
        sampleid='{sample}',
        random_nb=randint(1,1000),
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']


use rule purecn_cbs_Hclust as purecn_hmm_Hclust with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
        seg='results/cnvkit'+suffix+'/hmm/{sample}.seg'
    output: 'results/purecn'+suffix+'/hmm_Hclust/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/hmm_Hclust/{sample}.txt'
    log: 'logs/purecn'+suffix+'/hmm_Hclust/{sample}.log',
    params:
        cnvkit_method='hmm',
        purecn_method='Hclust',
        sampleid='{sample}',
        random_nb=randint(1,1000),
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']


use rule purecn_cbs_Hclust as purecn_hmm_PSCBS with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
        seg='results/cnvkit'+suffix+'/hmm/{sample}.seg',
        install='results/purecn/PSCBS_check'
    output: 'results/purecn'+suffix+'/hmm_PSCBS/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/hmm_PSCBS/{sample}.txt'
    log: 'logs/purecn'+suffix+'/hmm_PSCBS/{sample}.log',
    params:
        cnvkit_method='hmm',
        purecn_method='PSCBS',
        sampleid='{sample}',
        random_nb=randint(1,1000),
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']


use rule purecn_cbs_Hclust as purecn_cbs_none with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
        seg='results/cnvkit'+suffix+'/cbs/{sample}.seg'
    output: 'results/purecn'+suffix+'/cbs_none/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/cbs_none/{sample}.txt'
    log: 'logs/purecn'+suffix+'/cbs_none/{sample}.log',
    params:
        cnvkit_method='cbs',
        purecn_method='none',
        sampleid='{sample}',
        random_nb=randint(1,1000),
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']


use rule purecn_cbs_Hclust as purecn_hmm_none with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
        seg='results/cnvkit'+suffix+'/hmm/{sample}.seg'
    output: 'results/purecn'+suffix+'/hmm_none/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/hmm_none/{sample}.txt'
    log: 'logs/purecn'+suffix+'/hmm_none/{sample}.log',
    params:
        cnvkit_method='hmm',
        purecn_method='none',
        sampleid='{sample}',
        random_nb=randint(1,1000),
        suffix=suffix,
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        genome=config['genome_version']
