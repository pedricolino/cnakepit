rule download_ploidyNGS:
    output: 'workflow/ploidyNGS/ploidyNGS.py'
    log: 'logs/ploidyngs/download_tool.log'
    shell: 'git clone "https://github.com/diriano/ploidyNGS.git" "workflow/ploidyNGS" 2> {log}'

# requires modification, should be done with a single chromosome
rule ploidyNGS:
    input: 
        bam = BAMs_for_CNV_calling,
        tool = 'workflow/ploidyNGS/ploidyNGS.py'
    output: 'results/ploidyngs/{sample}/{sample}_MaxDepth100_MinCov0.tab'
    conda: env_prefix + 'cnv_calling' + env_suffix
    log: 'logs/ploidyngs/{sample}.log'
    resources: mem=lambda wildcards, attempt: '%dG' % (20 * attempt)
    shell: 
        '''
        cd workflow/ploidyNGS/ 
            ./ploidyNGS.py -b ../../{input.bam} -o {wildcards.sample} -g -c 100 > ../../{log}  2>&1   
            mkdir -p results/ploidyngs/{wildcards.sample} 
            mv {wildcards.sample}* results/ploidyngs/{wildcards.sample}/ 
            cd -
        '''

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

rule get_vep:
    output: config['vep']['data'] + '/homo_sapiens/109_GRCh37/'
    log: 'logs/get_vep.log'
    conda: env_prefix + 'vcf2maf' + env_suffix
    shell:
        '''
        mkdir -p {config['vep']['data']}
        vep_install -a cf -s homo_sapiens -y GRCh37 -c resources/vep/hg37 --CONVERT -n
        '''

rule filter_purecn_vcf:
    input: 'results/purecn'+suffix+'/cbs_Hclust/{sample}/{sample}.vcf'
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
        vep = config['vep']['data'] + '/homo_sapiens/109_GRCh37/'
    output:
        maf = 'results/vcf2maf/purecn'+suffix+'_cbs_none/{sample}.maf',
        # vcf = 'results/vcf2maf/purecn'+suffix+'_cbs_none/{sample}.vep.vcf'
    conda: env_prefix + 'vcf2maf' + env_suffix
    log: 'logs/vcf2maf/purecn'+suffix+'cbs_none/{sample}.log'
    params:
        sampleid = '{sample}',
        vep_path = config['vep']['executable'],
        vep_data = config['vep']['data'],
        folder = 'results/vcf2maf/purecn'+suffix+'_cbs_none',
    threads: 1
    resources: mem=lambda wildcards, attempt: '%dG' % (4 * 8 * attempt), # 4 GB per thread, 8 threads
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

rule merge_maf:
    input: expand('results/vcf2maf/purecn'+suffix+'_cbs_none/{sample}.maf', sample=samples_for_calling.index),
    output: 'results/vcf2maf/purecn'+suffix+'_cbs_none/allsamples.maf'
    shell:
        '''
        set -o pipefail
        cat {input} | egrep "^#|^Hugo_Symbol" | head -2 > {output} || true
        cat {input} | egrep -v "^#|^Hugo_Symbol" >> {output}
        '''

rule merge_seg:
    input: expand('results/purecn'+suffix+'/cbs_none/{sample}/{sample}_dnacopy.seg', sample=samples_for_calling.index),
    output: 'results/purecn'+suffix+'/cbs_none/allsamples.seg'
    shell:
        '''
        cat {input} | egrep "^ID" | head -1 > {output}
        cat {input} | egrep -v "^ID" >> {output}
        '''

rule gistic2:
    input:  'results/purecn'+suffix+'/cbs_none/allsamples.seg'
    output: 'results/gistic/purecn'+suffix+'_cbs_none/scores.gistic'
    conda: env_prefix + 'gistic' + env_suffix
    resources: mem=lambda wildcards, attempt: '%dG' % (8 * attempt)
    params:
        folder = 'results/gistic/purecn'+suffix+'_cbs_none/',
        refgene = config['gistic']['refgene'],
        conf = config['gistic']['conf_level'],
    benchmark: 'benchmarks/gistic/purecn'+suffix+'_cbs_none.txt'
    shell:
        '''
        rm -rf {params.folder}
        mkdir -p {params.folder}
        gistic2 \
            -seg {input} \
            -refgene {params.refgene} \
            -conf_level {params.conf} \
            -b {params.folder}
        '''