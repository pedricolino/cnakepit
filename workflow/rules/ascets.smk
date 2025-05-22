rule get_ascets:
    output:
        script = 'workflow/repos/ascets/run_ascets.R',
        ref = 'workflow/repos/ascets/genomic_arm_coordinates_hg19.txt' if config['genome_version'] == 'hg19' else 'workflow/repos/ascets/genomic_arm_coordinates_hg38.txt'
    params: config['ascets']['github_link']
    shell:
        '''
        git clone {params} workflow/repos/ascets
        '''

rule prepare_ascets_input:
    input:  'results/cnvkit'+suffix+'/cbs/{sample}.seg'
    output: 'results/ascets/purecn'+suffix+'_cbs/input/{sample}.seg'
    shell:
        '''
        cut -f 1-6 {input} > {output}
        '''

rule run_ascets:
    input:
        script = 'workflow/repos/ascets/run_ascets.R',
        cleaned_seg = 'results/ascets/purecn'+suffix+'_cbs/input/{sample}.seg',
        ref = 'workflow/repos/ascets/genomic_arm_coordinates_hg19.txt' if config['genome_version'] == 'hg19' else 'workflow/repos/ascets/genomic_arm_coordinates_hg38.txt'
    output:
        'results/ascets/purecn'+suffix+'_cbs/output/{sample}_arm_level_calls.txt'
    params:
        sampleid = 'results/ascets/purecn'+suffix+'_cbs/output/{sample}',
    conda: env_prefix + 'R' + env_suffix # needs tidyverse and data.table packages
    shell:
        '''
        Rscript {input.script} \
            -i {input.cleaned_seg} \
            -o {params.sampleid} \
            -c {input.ref}
        '''
