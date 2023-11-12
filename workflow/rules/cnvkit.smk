bedtargets = 'results/cnvkit/'+bedname+'_target.bed'
bedantitargets = 'results/cnvkit/'+bedname+'_antitarget.bed'
bedname=bedname

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_mappability:
    input:
        HTTP.remote("github.com/etal/cnvkit/raw/master/data/access-5k-mappable.hg19.bed", keep_local=True)
    output:
        config["mappability"]
    shell:
        "mv {input} {output}"

rule cnvkit_autobin:
    input:
        bams = expand("results/bam_sorted_bwa/{sample}_sorted.bam", sample=samples.index),
        targets = config["bed_w_chr"],
        access = config["mappability"],
    output:
        target = bedtargets,
        antitarget = bedantitargets,
    params:
        extra = '--method amplicon',
        samplenames = samples.index
    threads: 8
    log:
        "logs/cnvkit_autobin/log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/autobin',
    shell:
        'cnvkit.py autobin {input.bams} --targets {input.targets} --access {input.access} --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} {params.extra} 2> {log}'


rule cnvkit_coverage:
    input:
        bam = 'results/bam_sorted_bwa/{sample}_sorted.bam',
        targets = 'results/cnvkit/'+bedname+'_target.bed',
        antitargets = 'results/cnvkit/'+bedname+'_antitarget.bed',
    output:
        target_coverage = 'results/cnvkit/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/{sample}.antitargetcoverage.cnn',
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_coverage/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/coverage'
    shell:
        'cnvkit.py coverage {input.bam} {input.targets} -o {output.target_coverage} {params.extra} && '
        'cnvkit.py coverage {input.bam} {input.antitargets} -o {output.antitarget_coverage} {params.extra} 2> {log}'

rule cnvkit_ref_generic:
    input:
        fasta=config["ref"],
        targets = bedtargets
    output:
        FlatReference_cnn = 'results/cnvkit/FlatReference.cnn',
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_ref_generic/log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py reference -o {output.FlatReference_cnn} -f {input.fasta} -t {input.targets} {params.extra} 2> {log}'

rule cnvkit_fix:
    input:
        target_coverage = 'results/cnvkit/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/{sample}.antitargetcoverage.cnn',
        reference = 'results/cnvkit/FlatReference.cnn',
    output:
        'results/cnvkit/{sample}.cnr'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_fix/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/fix'
    shell:
        'cnvkit.py fix {input.target_coverage} {input.antitarget_coverage} {input.reference} -o {output} {params.extra} 2> {log}'

rule cnvkit_segment:
    input:
        copy_ratios = 'results/cnvkit/{sample}.cnr',
    output:
        'results/cnvkit/{sample}.cns',
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_segment/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/segment'
    shell:
        'cnvkit.py segment {input.copy_ratios} -o {output} {params.extra} 2> {log}'

rule cnvkit_scatter:
    input:
        copy_ratio = 'results/cnvkit/{sample}.cnr',
        segment = 'results/cnvkit/{sample}.cns',
    output:
        'results/cnvkit/{sample}_scatter.cnv.pdf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Plot segment lines in this color. value can be any string
        # accepted by matplotlib, e.g. 'red' or '#CC0000'
        #segment_color = '',
        # Plot title. [Default: sample ID, from filename or -i]
        #title = '',
    threads: 8
    log:
        "logs/cnvkit_scatter/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/scatter'
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram:
    input:
        copy_ratio = 'results/cnvkit/{sample}.cnr',
        segment = 'results/cnvkit/{sample}.cns',
    output:
        'results/cnvkit/{sample}_diagram.cnv.pdf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Plot segment lines in this color. value can be any string
        # accepted by matplotlib, e.g. 'red' or '#CC0000'
        #segment_color = '',
        # Plot title. [Default: sample ID, from filename or -i]
        #title = '',
    threads: 8
    log:
        "logs/cnvkit_diagram/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/diagram'
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap:
    input:
        segment = 'results/cnvkit/{sample}.cns'
    output:
        'results/cnvkit/{sample}_heatmap.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_heatmap/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segment} {params.extra} 2> {log}'

rule export_seg:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit/{sample}.cns',
    output:
        'results/cnvkit/{sample}.seg'
    params:
        extra = '--enumerate-chroms',
    threads: 8
    log:
        "logs/export_seg/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'


# rule cnvkit-batch-amplicon-nocontrol:
#     input:
#         fasta=config["ref"],
#         map="results/bam_sorted_bwa/{sample}_sorted.bam",
#         access=config["mappability"]
#     output:
#         targetbed="results/cnvkit/{bedname}.target.bed",
#         antitargetbed="results/cnvkit/{bedname}.antitarget.bed",
#         scatter="{sample}_sorted-scatter.png",
#         diagram="{sample}_sorted-diagram.pdf",
#         targetcoverage="{sample}_sorted.targetcoverage.cnn",
#         cns="{sample}_sorted.cns",
#         cnr="{sample}_sorted.cnr",
#         cnr="{sample}_sorted.cnr",
#         callcns="{sample}_sorted.call.cns",
#         bintestcns="{sample}_sorted.bintest.cns",
#         antitargetcoveragecnn="{sample}_sorted.antitargetcoverage.cnn",
#         callcns="{sample}_sorted.call.cns"
#     message:
#         "Running cnvkit batch amplicon no control for {sample}"
#     threads: 8
#     log:
#         "logs/cnvkit/{sample}_sorted.log"
#     conda:
#         "envs/primary_env.yaml"
#     shell:
#         ""
        

############# added rules for additional segmentation method hmm

rule cnvkit_segment_hmm:
    input:
        copy_ratios = 'results/cnvkit/{sample}.cnr',
    output:
        'results/cnvkit_hmm/{sample}.cns',
    params:
        extra = '-m hmm',
    threads: 8
    log:
        "logs/cnvkit_segment_hmm/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/segment'
    shell:
        'cnvkit.py segment {input.copy_ratios} -o {output} {params.extra} 2> {log}'


rule cnvkit_scatter_hmm:
    input:
        copy_ratio = 'results/cnvkit/{sample}.cnr',
        segment = 'results/cnvkit_hmm/{sample}.cns',
    output:
        'results/cnvkit_hmm/{sample}_scatter.cnv.pdf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Plot segment lines in this color. value can be any string
        # accepted by matplotlib, e.g. 'red' or '#CC0000'
        #segment_color = '',
        # Plot title. [Default: sample ID, from filename or -i]
        #title = '',
    threads: 8
    log:
        "logs/cnvkit_scatter_hmm/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/scatter'
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram_hmm:
    input:
        copy_ratio = 'results/cnvkit/{sample}.cnr',
        segment = 'results/cnvkit_hmm/{sample}.cns',
    output:
        'results/cnvkit_hmm/{sample}_diagram.cnv.pdf'
    params:
        # Optional parameters. Omit if unneeded.
        extra = '',
        # Plot segment lines in this color. value can be any string
        # accepted by matplotlib, e.g. 'red' or '#CC0000'
        #segment_color = '',
        # Plot title. [Default: sample ID, from filename or -i]
        #title = '',
    threads: 8
    log:
        "logs/cnvkit_diagram_hmm/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    # wrapper:
    #     'http://dohlee-bio.info:9193/cnvkit/diagram'
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap_hmm:
    input:
        segment = 'results/cnvkit_hmm/{sample}.cns'
    output:
        'results/cnvkit_hmm/{sample}_heatmap.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit_heatmap_hmm/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segment} {params.extra} 2> {log}'

rule export_seg_hmm:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit_hmm/{sample}.cns',
    output:
        'results/cnvkit_hmm/{sample}.seg'
    params:
        extra = '--enumerate-chroms',
    threads: 8
    log:
        "logs/export_seg_hmm/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'
