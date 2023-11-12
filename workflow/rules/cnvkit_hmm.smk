rule cnvkit_segment_hmm:
    input:
        copy_ratios = 'results/cnvkit/general/{sample}.cnr',
    output:
        'results/cnvkit/hmm/{sample}.cns',
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/hmm/segment/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py segment {input.copy_ratios} -o {output} {params.extra} 2> {log}'

rule cnvkit_scatter_hmm:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
    output:
        'results/cnvkit/hmm/{sample}_scatter.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/hmm/scatter/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram_hmm:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
    output:
        'results/cnvkit/hmm/{sample}_diagram.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/hmm/diagram/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap_hmm:
    input:
        segment = 'results/cnvkit/hmm/{sample}.cns'
    output:
        'results/cnvkit/hmm/{sample}_heatmap.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/hmm/heatmap/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segment} {params.extra} 2> {log}'

rule export_seg_hmm:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit/hmm/{sample}.cns',
    output:
        'results/cnvkit/hmm/{sample}.seg'
    params:
        extra = '--enumerate-chroms',
    threads: 8
    log:
        "logs/cnvkit/hmm/export_seg/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'
