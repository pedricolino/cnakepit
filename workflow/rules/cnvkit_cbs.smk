rule cnvkit_segment_cbs:
    input:
        copy_ratios = 'results/cnvkit/general/{sample}.cnr',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/cbs/{sample}.cns',
    benchmark:
        'benchmarks/cnvkit/cbs_segment/{sample}.txt'
    params:
        extra = '-m cbs --drop-low-coverage --smooth-cbs'
    threads: 8
    log:
        "logs/cnvkit/cbs/segment/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py segment {input.copy_ratios} --vcf {input.germline_vcf} -o {output} --processes {threads} {params.extra} 2> {log}'

rule cnvkit_scatter_cbs:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/cbs/{sample}.cns',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/cbs/{sample}_scatter.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/cbs/{sample}_scatter.txt'
    params:
        extra = '',
    log:
        "logs/cnvkit/cbs/scatter/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram_cbs:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/cbs/{sample}.cns',
    output:
        'results/cnvkit/cbs/{sample}_diagram.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/cbs/{sample}_diagram.txt'
    params:
        extra = '',
    log:
        "logs/cnvkit/cbs/diagram/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap_cbs:
    input:
        segments = expand("results/cnvkit/cbs/{sample}.cns", sample=samples.index)
    output:
        'results/cnvkit/cbs/heatmap.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/cbs/heatmap.txt'
    log:
        "logs/cnvkit/cbs/heatmap.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segments} &> {log}'

rule export_seg_cbs:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit/cbs/{sample}.cns',
    output:
        'results/cnvkit/cbs/{sample}.seg'
    benchmark:
        'benchmarks/cnvkit/cbs/export_seg/{sample}.txt'
    params:
        extra = '--enumerate-chroms',
    log:
        "logs/cnvkit/cbs/export_seg/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'

rule cnvkit_call_cbs:
    input:
        cns = 'results/cnvkit/cbs/{sample}.cns',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/cbs/{sample}.icns' # integer copy number
    benchmark:
        'benchmarks/cnvkit/cbs/{sample}.call.txt'
    params:
        extra = '--drop-low-coverage', # https://cnvkit.readthedocs.io/en/stable/tumor.html
        # maybe include --center median
    log:
        "logs/cnvkit/cbs/call/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py call {input.cns} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'