rule cnvkit_segment_cbs:
    input:
        copy_ratios = 'results/cnvkit/general/{sample}.cnr',
    output:
        'results/cnvkit/cbs/{sample}.cns',
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/cbs/segment/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py segment {input.copy_ratios} -o {output} {params.extra} 2> {log}'

rule cnvkit_scatter_cbs:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/cbs/{sample}.cns',
    output:
        'results/cnvkit/cbs/{sample}_scatter.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/cbs/scatter/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram_cbs:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/cbs/{sample}.cns',
    output:
        'results/cnvkit/cbs/{sample}_diagram.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/cbs/diagram/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap_cbs:
    input:
        segment = 'results/cnvkit/cbs/{sample}.cns'
    output:
        'results/cnvkit/cbs/{sample}_heatmap.cnv.pdf'
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/cbs/heatmap/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segment} {params.extra} 2> {log}'

rule export_seg_cbs:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit/cbs/{sample}.cns',
    output:
        'results/cnvkit/cbs/{sample}.seg'
    params:
        extra = '--enumerate-chroms',
    threads: 8
    log:
        "logs/cnvkit/cbs/export_seg/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'
