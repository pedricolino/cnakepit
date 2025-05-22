# merged file required by gistic2
rule merge_seg:
    input: expand('results/purecn'+suffix+'/cbs_none/{sample}/{sample}_dnacopy.seg', sample=samples_for_calling.index),
    output: 'results/purecn'+suffix+'/cbs_none/allsamples.seg'
    shell:
        '''
        cat {input} | egrep "^ID" | head -1 > {output}
        cat {input} | egrep -v "^ID" >> {output}
        '''

rule gistic2:
    input:  'results/purecn'+suffix+'/cbs_none/allsamples.seg',
    output:
        'results/gistic/purecn'+suffix+'_cbs_none/scores.gistic',
        'results/gistic/purecn'+suffix+'_cbs_none/all_thresholded.by_genes.txt'
    conda: env_prefix + 'gistic' + env_suffix
    resources: mem=lambda wildcards, attempt: '%dG' % (8 * attempt)
    params:
        folder = 'results/gistic/purecn'+suffix+'_cbs_none/',
        refgene = config['reference'+'_'+config['genome_version']]["gistic_refgene"],
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
            -b {params.folder} \
            -savegene 1 \
            -save_data_files 0 \
            -v 30
        '''
		

# create the 3 files required by cbioportal (at least CUBI's version)

rule gistic2_results_to_cbioportal_format:
    input:  'results/gistic/purecn'+suffix+'_cbs_none/all_thresholded.by_genes.txt'
    output: 'results/purecn'+suffix+'/cbs_none/for_cbioportal/data_cna_gistic.txt'
    shell:
        '''
        (
        head -n 1 {input} | awk 'BEGIN{{OFS="\t"}} {{$1="Hugo_Symbol"; $2="Entrez_Gene_Id"; print}}'
        tail -n +2 {input} | sed -E $'s/\\|chr[0-9]+//' | awk 'BEGIN{{OFS="\t"}} {{for(i=1;i<=NF;i++) if(i!=3) printf "%s%s", $i, (i==NF||i==NF-1&&3==NF ? RS : OFS)}}'
        ) > {output}
        '''

rule purecn_results_to_cbioportal_log2_format:
    input:  expand('results/purecn'+suffix+'/cbs_none/{sample}/{sample}_genes.csv', sample=samples_for_calling.index)
    output: 'results/purecn'+suffix+'/cbs_none/for_cbioportal/data_cna_log2.txt'
    conda: env_prefix + 'R' + env_suffix
    shell:
        '''
        Rscript workflow/scripts/purecn_to_cbioportal_log2_format.R {output} {input}
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