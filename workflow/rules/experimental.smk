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