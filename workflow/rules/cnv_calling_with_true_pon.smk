# These rules require pre-made pon.rds and interval.rds files
# The original commands from the bash scripts are kept in the comments

# Rscript $PURECN/Coverage.R \
#     --out-dir $OUT/$SAMPLEID \
#     --bam $BAMFOLDER/$SAMPLEID.bam \
#     --intervals $INTERVALS.no_chr.txt \
#     --force

rule coverage:
    input:
        bam = BAMs_for_CNV_calling,
        intervals = config['pon_rds']['intervals_file']
    output: 'results/cnv_calling_with_true_pon/coverage/{sample}/{sample}_coverage_loess.txt.gz'
    params:
        sampleid='{sample}',
    conda: env_prefix + 'cnv_calling' + env_suffix
    log: 'logs/cnv_calling_with_true_pon/coverage/{sample}.log'
    benchmark: 'benchmarks/cnv_calling_with_true_pon/coverage/{sample}.tsv'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
    threads: 1
    shell:
        '''
        PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
        Rscript $PURECN/Coverage.R \
        --out-dir results/cnv_calling_with_true_pon/coverage/{params.sampleid} \
        --bam {input.bam} \
        --intervals {input.intervals} 2> {log}
        '''

# zcat $OUT/$SAMPLEID/${SAMPLEID}_coverage_loess.txt.gz |
#     awk 'NR>1 {print "chr"$0; next} 1' >$OUT/$SAMPLEID/${SAMPLEID}_coverage_loess_w_chr.txt.gz

if config['pon_rds']['different_contigs']:
    # add 'chr' to the contig lines to make it compatible with the used reference genome
    rule add_chr_to_coverage:
        input:
            'results/cnv_calling_with_true_pon/coverage/{sample}/{sample}_coverage_loess.txt.gz'
        output:
            'results/cnv_calling_with_true_pon/coverage/{sample}/{sample}_coverage_loess_w_chr.txt.gz'
        shell:
            """
            zcat {input} | awk 'NR>1 {{print "chr"$0; next}} 1' >{output}
            """

    # gatk --java-options "-Xmx80G" LiftoverVcf \
    #     -I $MUTECTFOLDER/$SAMPLEID.vcf.gz \
    #     -O $OUT/$SAMPLEID/${SAMPLEID}_hg19.vcf.gz \
    #     -C b37tohg19.chain \
    #     --REJECT $OUT/$SAMPLEID/${SAMPLEID}_rejected_variants.vcf.gz \
    #     --REFERENCE_SEQUENCE $FASTA \
    #     --TMP_DIR ~/scratch/tmp

    # liftover the vcf file from b37 (required by Mutect2) to hg19 (required by custom pon.rds)
    rule liftover_vcf:
        input:
            vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
            chain = config['pon_rds']['chain_file'],
            fasta = config['pon_rds']['fasta']
        output:
            lifted_over='results/cnv_calling_with_true_pon/vcf/{sample}/{sample}_hg19.vcf.gz',
            rejected='results/cnv_calling_with_true_pon/vcf/{sample}/{sample}_rejected_variants.vcf.gz'
        params:
            sampleid='{sample}',
        resources:
            mem=lambda wildcards, attempt: '%dG' % (8 * attempt),
        conda: env_prefix + 'cnv_calling' + env_suffix
        log: 'logs/cnv_calling_with_true_pon/liftover_vcf/{sample}.log'
        benchmark: 'benchmarks/cnv_calling_with_true_pon/liftover_vcf/{sample}.tsv'
        threads: 1
        shell:
            "gatk --java-options -Xmx{resources.mem} LiftoverVcf "
                "-I {input.vcf_filt} "
                "-O {output.lifted_over} "
                "-C {input.chain} "
                "--REJECT {output.rejected} "
                "--REFERENCE_SEQUENCE {input.fasta} "
                "--TMP_DIR ~/scratch/tmp 2> {log}"


rule download_sv_blacklist:
    input:
        HTTP.remote(config['sv_blacklist'+'_'+config['genome_version']]['link'], keep_local=True)
    output:
        config['sv_blacklist'+'_'+config['genome_version']]['bed']
    benchmark:
        'benchmarks/cnvkit/general/download_sv_blacklist.txt'
    log:
        'logs/cnvkit/general/download_sv_blacklist/log',
    shell:
        'mv {input} {output}'


# Rscript $PURECN/PureCN.R --out $OUT/$SAMPLEID \
#     --tumor $OUT/$SAMPLEID/${SAMPLEID}_coverage_loess_w_chr.txt.gz \
#     --sampleid $SAMPLEID \
#     --vcf $OUT/$SAMPLEID/${SAMPLEID}_hg19.vcf.gz \
#     --stats-file $MUTECTFOLDER/$SAMPLEID.vcf.gz.stats \
#     --fun-segmentation PSCBS \
#     --normaldb $WHICHREF/normaldb.rds \
#     --mapping-bias-file $WHICHREF/mapping-bias-file.rds \
#     --intervals $INTERVALS \
#     --snp-blacklist $BLACKLIST \
#     --genome hg19 \
#     --force --post-optimize --seed 123

# parameters as they are in the pipeline from patho
rule purecn_true_pon_cbs:
    input:
        coverage = 'results/cnv_calling_with_true_pon/coverage/{sample}/{sample}_coverage_loess_w_chr.txt.gz' if config['pon_rds']['different_contigs'] else 'results/cnv_calling_with_true_pon/coverage/{sample}/{sample}_coverage_loess.txt.gz',
        vcf = 'results/cnv_calling_with_true_pon/vcf/{sample}/{sample}_hg19.vcf.gz' if config['pon_rds']['different_contigs'] else 'results/mutect2/filtered/{sample}_filtered.vcf.gz',
        stats = 'results/mutect2/unfiltered/{sample}.vcf.gz.stats',
        normaldb = config['pon_rds']['rds'],
        mapping_bias = config['pon_rds']['mapping_bias'] if config['pon_rds']['mapping_bias_available'] else '',
        intervals = config['pon_rds']['intervals_file'],
        blacklist = config['sv_blacklist'+'_'+config['genome_version']]['bed'],
    # output: 'results/cnv_calling_with_true_pon/purecn/CBS/{sample}/{sample}.csv'
    output: 'results/cnv_calling_with_true_pon/purecn/CBS/{sample}/{sample}.csv'
    params:
        sampleid='{sample}',
        genome=config['genome_version'],
        minpurity=0.1,
        maxploidy=15,
        maxcopynumber=20,
        mincoverage=50,
        random_nb=randint(1,1000), # uses a different seed on retry in case of failure
        sex= '' if not config['sex']['hard_code'] else '--sex '+config['sex']['sex'],
        bootstrapping=100,
        purecn_method='CBS',
        mappingbias= '' if not config['pon_rds']['mapping_bias_available'] else '--mapping-bias-file',
    resources: mem=lambda wildcards, attempt: '%dG' % (4 * 8 * attempt), # 4 GB per thread, 8 threads
    threads: 8   
    conda: env_prefix + 'cnv_calling' + env_suffix
    # log: 'logs/cnv_calling_with_true_pon/purecn/cbs/{sample}.log'
    log: 'logs/cnv_calling_with_true_pon/purecn/CBS/{sample}.log'
    # benchmark: 'benchmarks/cnv_calling_with_true_pon/purecn/cbs/{sample}.tsv'
    benchmark: 'benchmarks/cnv_calling_with_true_pon/purecn/CBS/{sample}.tsv'
    shell:
        """
        PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')
        mkdir -p results/cnv_calling_with_true_pon/purecn/{params.purecn_method}/{params.sampleid}
        Rscript $PURECN/PureCN.R \
            --tumor {input.coverage} \
            --sampleid {params.sampleid} \
            --vcf {input.vcf} \
            --stats-file {input.stats} \
            --fun-segmentation {params.purecn_method} \
            --normaldb {input.normaldb} \
            {params.mappingbias} {input.mapping_bias} \
            --intervals {input.intervals} \
            --snp-blacklist {input.blacklist} \
            --genome {params.genome} \
            --out results/cnv_calling_with_true_pon/purecn/{params.purecn_method}/{params.sampleid}/{params.sampleid} \
            --min-purity {params.minpurity} \
            --max-ploidy {params.maxploidy} \
            --max-copy-number {params.maxcopynumber} \
            --min-total-counts {params.mincoverage} \
            --cores {threads} \
            --seed {params.random_nb} \
            --bootstrap-n {params.bootstrapping} \
            {params.sex} \
            --force --post-optimize --seed 123 2> {log}
            """
