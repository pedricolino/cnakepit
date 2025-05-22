#-------------- Preparations for vcf2maf

rule download_vcf2maf:
    output: 'resources/mskcc-vcf2maf-f6d0c40/vcf2maf.pl'
    log: 'logs/vcf2maf/download_tool.log'
    shell:
        '''
        curl -L -o resources/mskcc-vcf2maf.tar.gz https://api.github.com/repos/mskcc/vcf2maf/tarball/v1.6.22 
        cd resources/
        tar -zxf mskcc-vcf2maf.tar.gz
        rmz mskcc-vcf2maf.tar.gz
        '''

# Hg19 is hardcoded in this and the vcf2maf rule. Fix this later.
rule get_vep:
    output: config['reference'+'_'+config['genome_version']]["vep_data"] + '/homo_sapiens/109_GRCh37/'
    log: 'logs/get_vep.log'
    conda: env_prefix + 'vcf2maf' + env_suffix
    params: config['reference'+'_'+config['genome_version']]["vep_data"]
    shell:
        '''
        mkdir -p {params}
        vep_install -a cf -s homo_sapiens -y GRCh37 -c resources/vep/hg37 --CONVERT -n
        '''

#-------------- To analyze only somatic variants as classified by mutect2

rule filter_purecn_vcf:
    input: 'results/purecn'+suffix+'/cbs_none/{sample}/{sample}.vcf'
    output: 'results/purecn'+suffix+'/cbs_none/{sample}/{sample}.filtered.vcf'
    conda: env_prefix + 'cnv_calling' + env_suffix
    shell:
        '''
        bcftools view -f PASS {input} > {output}
        '''

rule vcf2maf:
    input:  
        vcf = 'results/purecn'+suffix+'/cbs_none/{sample}/{sample}.filtered.vcf',
        script = 'resources/mskcc-vcf2maf-f6d0c40/vcf2maf.pl',
        ref = ref_file,
        vep = config['reference'+'_'+config['genome_version']]["vep_data"] + '/homo_sapiens/109_GRCh37/'
    output:
        maf = 'results/vcf2maf/purecn'+suffix+'_cbs_none/{sample}.maf',
        # vcf = 'results/vcf2maf/purecn'+suffix+'_cbs_none/{sample}.filtered.vep.vcf' # not necessary
    conda: env_prefix + 'vcf2maf' + env_suffix
    log: 'logs/vcf2maf/purecn'+suffix+'cbs_none/{sample}.log'
    params:
        sampleid = '{sample}',
        vep_path = config['vep']['executable'],
        vep_data = config['reference'+'_'+config['genome_version']]["vep_data"],
        folder = 'results/vcf2maf/purecn'+suffix+'_cbs_none',
    threads: 1
    resources: mem=lambda wildcards, attempt: '%dG' % (4 * 8 * attempt)
    shell:
        '''
        perl {input.script} \
            --input-vcf {input.vcf} \
            --output-maf {output.maf} \
            --tmp-dir {params.folder} \
            --ref-fasta {input.ref} \
            --vep-path {params.vep_path} \
            --vep-data {params.vep_data} \
            --tumor-id {params.sampleid} \
            --vep-overwrite \
            --vep-forks {threads} 2> {log}
        '''

#-------------- To analyze only somatic variants as classified by mutect2

# TBD

#-------------- To analyze unfiltered (germline+somatic+errors) variants from mutect2

# vcf2maf requires uncompressed vcf format
rule uncompress_vcf:
    input: 'results/mutect2/filtered/{sample}_filtered.vcf.gz'
    output: temp('results/mutect2/filtered/{sample}_filtered.vcf') # delete when not needed anymore
    shell:
        '''
        gunzip -c {input} > {output}
        '''

use rule vcf2maf as vcf2maf_unfiltered with:
    input:
        vcf = 'results/mutect2/filtered/{sample}_filtered.vcf',
        script = 'resources/mskcc-vcf2maf-f6d0c40/vcf2maf.pl',
        ref = ref_file,
        vep = config['reference'+'_'+config['genome_version']]["vep_data"] + '/homo_sapiens/109_GRCh37/'
    output: maf = 'results/vcf2maf/unfiltered/{sample}_unfiltered.maf'
    log: 'logs/vcf2maf/unfiltered/{sample}.log'
    params:
        sampleid = '{sample}',
        vep_path = config['vep']['executable'],
        vep_data = config['reference'+'_'+config['genome_version']]["vep_data"],
        folder = 'results/vcf2maf/unfiltered',

use rule merge_maf as merge_maf_unfiltered with:
    input: expand('results/vcf2maf/unfiltered/{sample}_unfiltered.maf', sample=samples_for_calling.index),
    output: 'results/vcf2maf/unfiltered/allsamples_unfiltered.maf'
