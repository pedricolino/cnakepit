import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# produced bed files
my_targets = 'results/cnvkit/general/my_targets.bed'
my_antitargets = 'results/cnvkit/general/my_antitargets.bed'

autobin_targets = 'results/cnvkit/general/'+bedname+'_target.bed'
autobin_antitargets = 'results/cnvkit/general/'+bedname+'_antitarget.bed'

bedname=bedname


# if second round, create 2 long strings with all PON (anti-/)target coverage files
if config['pon']['second_run_with_pon'] == True:

    # Open the file and read the lines
    with open(config['pon']['reference_sample_names'], 'r') as file:
        names = file.readlines()

    # Remove newline characters

    names = [name.strip() for name in names]

    # Format the sample names
    targets_str = ['results/cnvkit/general/'+name+'.targetcoverage.cnn' for name in names]
    antitargets_str = ['results/cnvkit/general/'+name+'.antitargetcoverage.cnn' for name in names]

    # Join the formatted strings into a single string with spaces
    long_str_target = ' '.join(targets_str)
    long_str_antitarget = ' '.join(antitargets_str)

    #   create reference track from the PON
    rule cnvkit_create_panel_of_normals:
        input:
            fasta=config['reference']['fasta'],
            all_target_cnns = expand('results/cnvkit/general/{sample}.targetcoverage.cnn', sample=samples.index),
            all_antitarget_cnns = expand('results/cnvkit/general/{sample}.antitargetcoverage.cnn', sample=samples.index)
        output:
            pon = 'results/cnvkit'+suffix+'/general/pon.cnn',
        benchmark: 'benchmarks/cnvkit'+suffix+'/general/cnvkit_create_panel_of_normals.txt'
        params:
            extra = '-c',
            sex= "" if not config['sex']['hard_code'] else '--sample-sex '+config['sex']['sex'].lower(),
            amplicon = '--no-edge' if config['amplicon'] else '',
            target_cnn = long_str_target,
            antitarget_cnn = long_str_antitarget,
        resources:
            mem=lambda wildcards, attempt: '%dG' % (16 * attempt),
        log:
            'logs/cnvkit'+suffix+'/general/cnvkit_create_panel_of_normals/log',
        conda:
            '../envs/cnv_calling.yaml'
        shell:
            'cnvkit.py reference -o {output.pon} -f {input.fasta} {params.extra} {params.sex} {params.amplicon} {params.target_cnn} {params.antitarget_cnn} 2> {log}'



# use remote file instead of cnvkit.py access output which causes problems
rule download_mappability:
    input:
        HTTP.remote(config['mappability']['link'], keep_local=True)
    output:
        config['mappability']['bed']
    benchmark:
        'benchmarks/cnvkit/general/download_mappability.txt'
    log:
        'logs/cnvkit/general/download_mappability/log',
    shell:
        'mv {input} {output}'

rule download_sv_blacklist:
    input:
        HTTP.remote(config['sv_blacklist']['link'], keep_local=True)
    output:
        config['sv_blacklist']['bed']
    benchmark:
        'benchmarks/cnvkit/general/download_sv_blacklist.txt'
    log:
        'logs/cnvkit/general/download_sv_blacklist/log',
    shell:
        'mv {input} {output}'

rule cnvkit_access:
    input:
        ref = config['reference']['fasta'],
        sv_blacklist = config['sv_blacklist']['bed'],
    output:
        'results/cnvkit/general/access.bed'
    benchmark:
        'benchmarks/cnvkit/general/access.txt'
    log:
        'logs/cnvkit/general/access/log',
    params:
        extra = '',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py access {input.ref} --exclude {input.sv_blacklist} -o {output} {params.extra} 2> {log}'

rule cnvkit_target:
    input:
        config['panel_design'],
    output:
        my_targets,
    benchmark:
        'benchmarks/cnvkit/general/target.txt'
    log:
        'logs/cnvkit/general/target/log',
    params:
        extra = '--split',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py target {input} -o {output} {params.extra} 2> {log}'

rule cnvkit_antitarget:
    input:
        bed = config['panel_design'],
        access = mappability
    output:
        my_antitargets,
    benchmark:
        'benchmarks/cnvkit/general/antitarget.txt'
    log:
        'logs/cnvkit/general/antitarget/log',
    params:
        extra = '',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py antitarget {input.bed} --access {input.access} -o {output} {params.extra} 2> {log}'

rule cnvkit_autobin:
    input:
        bams = expand(BAMs_for_CNV_calling, sample=samples.index), # maybe use only fraction of samples here?
        targets = my_targets,
        access = mappability,
        antitarget = my_antitargets,
    output:
        target = autobin_targets,
        antitarget = autobin_antitargets,
    resources:
        mem_mb=8000 # otherwise 4TB for 600 samples with default Snakemake...
    benchmark:
        'benchmarks/cnvkit/general/autobin.txt'
    params:
        extra = '',
        method = 'amplicon' if config['amplicon'] else 'hybrid',
        samplenames = samples.index
    log:
        'logs/cnvkit/general/autobin/log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py autobin {input.bams} --targets {input.targets} --access {input.access} --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} --method {params.method} {params.extra} 2> {log}'


rule cnvkit_coverage:
    input:
        bam = BAMs_for_CNV_calling,
        targets = autobin_targets,
        antitargets = autobin_antitargets,
    output:
        target_coverage = 'results/cnvkit/general/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/general/{sample}.antitargetcoverage.cnn',
    benchmark: 'benchmarks/cnvkit/general/coverage/{sample}.txt'
    params:
        extra = '',
    threads: 8
    log:
        'logs/cnvkit/general/coverage/{sample}.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py coverage {input.bam} {input.targets} --processes {threads} -o {output.target_coverage} {params.extra} && '
        'cnvkit.py coverage {input.bam} {input.antitargets} --processes {threads} -o {output.antitarget_coverage} {params.extra} 2> {log}'

rule cnvkit_generic_ref:
    input:
        fasta=config['reference']['fasta'],
        targets = autobin_targets,
        antitargets = autobin_antitargets,
    output:
        FlatReference_cnn = 'results/cnvkit/general/FlatReference.cnn',
    benchmark: 'benchmarks/cnvkit/general/ref_generic.txt'
    params:
        extra = '',
        sex= "" if not config['sex']['hard_code'] else '--sample-sex '+config['sex']['sex'].lower(),
        amplicon = '--no-edge' if config['amplicon'] else '',
    log:
        'logs/cnvkit/general/ref_generic/log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py reference -o {output.FlatReference_cnn} -f {input.fasta} -t {input.targets} -a {input.antitargets} {params.extra} {params.sex} {params.amplicon} 2> {log}'

rule cnvkit_fix:
    input:
        target_coverage = 'results/cnvkit/general/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/general/{sample}.antitargetcoverage.cnn',
        reference = 'results/cnvkit/general/FlatReference.cnn' if not config['pon']['second_run_with_pon'] else 'results/cnvkit'+suffix+'/general/pon.cnn',
    params:
        extra = '',
        amplicon = '--no-edge' if config['amplicon'] else '',
    output: 'results/cnvkit'+suffix+'/general/{sample}.cnr'
    benchmark: 'benchmarks/cnvkit'+suffix+'/general/fix/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/general/fix/{sample}.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py fix {input.target_coverage} {input.antitarget_coverage} {input.reference} -o {output} {params.extra} {params.amplicon} 2> {log}'


# CNVkit segmentation and downstream rules

rule cnvkit_segment_cbs:
    input:
        copy_ratios = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs_segment/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/segment/{sample}.log',
    params: extra = '-m cbs --drop-low-coverage --smooth-cbs'
    threads: 8
    resources: mem=lambda wildcards, attempt: '%dG' % (4 * attempt),
    conda: '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py segment {input.copy_ratios} --vcf {input.germline_vcf} -o {output} --processes {threads} {params.extra} 2> {log}'


rule cnvkit_scatter_cbs:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}_scatter.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs/{sample}_scatter.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/scatter/{sample}.log',
    params: extra = '',
    conda: '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'


rule cnvkit_diagram_cbs:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
    output:
        'results/cnvkit'+suffix+'/cbs/{sample}_diagram.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit'+suffix+'/cbs/{sample}_diagram.txt'
    params:
        extra = '',
    log:
        'logs/cnvkit'+suffix+'/cbs/diagram/{sample}.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'


rule cnvkit_heatmap_cbs:
    input:
        segments = expand('results/cnvkit'+suffix+'/cbs/{sample}.cns', sample=samples.index)
    output:
        'results/cnvkit'+suffix+'/cbs/heatmap.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit'+suffix+'/cbs/heatmap.txt'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        slurm_partition = 'medium'
    log:
        'logs/cnvkit'+suffix+'/cbs/heatmap.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py heatmap -o {output} {input.segments} &> {log}'


rule export_seg_cbs:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
    output:
        'results/cnvkit'+suffix+'/cbs/{sample}.seg'
    benchmark:
        'benchmarks/cnvkit'+suffix+'/cbs/export_seg/{sample}.txt'
    params:
        extra = '--enumerate-chroms',
    log:
        'logs/cnvkit'+suffix+'/cbs/export_seg/{sample}.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'


rule cnvkit_call_cbs:
    input:
        cns = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output:
        'results/cnvkit'+suffix+'/cbs/{sample}.icns' # integer copy number
    benchmark:
        'benchmarks/cnvkit'+suffix+'/cbs/{sample}.call.txt'
    params:
        extra = '--drop-low-coverage', # https://cnvkit.readthedocs.io/en/stable/tumor.html
        # maybe include --center median
    log:
        'logs/cnvkit'+suffix+'/cbs/call/{sample}.log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py call {input.cns} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'


# Do the same with the HMM segmentation algorithm

use rule cnvkit_segment_cbs as cnvkit_segment_hmm with:
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm_segment/{sample}.txt'
    params: extra = '-m hmm --drop-low-coverage',
    log: 'logs/cnvkit'+suffix+'/hmm/segment/{sample}.log',

use rule cnvkit_scatter_cbs as cnvkit_scatter_hmm with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/hmm/{sample}_scatter.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}_scatter.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/scatter/{sample}.log',

use rule cnvkit_diagram_cbs as cnvkit_diagram_hmm with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
    output: 'results/cnvkit'+suffix+'/hmm/{sample}_diagram.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}_diagram.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/diagram/{sample}.log',

use rule cnvkit_heatmap_cbs as cnvkit_heatmap_hmm with:
    input: segments = expand('results/cnvkit'+suffix+'/hmm/{sample}.cns', sample=samples.index)
    output: 'results/cnvkit'+suffix+'/hmm/heatmap.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/heatmap.txt'

use rule export_seg_cbs as export_seg_hmm with:
    input: cns = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.seg'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/export_seg/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/export_seg/{sample}.log',

use rule cnvkit_call_cbs as cnvkit_call_hmm with:
    input:
        cns = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.icns' # integer copy number
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}.call.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/call/{sample}.log',