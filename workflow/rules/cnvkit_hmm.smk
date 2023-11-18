rule cnvkit_segment_hmm:
    input:
        copy_ratios = 'results/cnvkit/general/{sample}.cnr',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/hmm/{sample}.cns',
    benchmark:
        'benchmarks/cnvkit/hmm_segment/{sample}.txt'
    params:
        extra = '-m hmm',
    threads: 8
    log:
        "logs/cnvkit/hmm/segment/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py segment {input.copy_ratios} --vcf {input.germline_vcf} -o {output} --processes {threads} {params.extra} 2> {log}'

rule cnvkit_scatter_hmm:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/hmm/{sample}_scatter.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/hmm/{sample}_scatter.txt'
    params:
        extra = '',
    log:
        "logs/cnvkit/hmm/scatter/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py scatter {input.copy_ratio} --segment {input.segment} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'

rule cnvkit_diagram_hmm:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
    output:
        'results/cnvkit/hmm/{sample}_diagram.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/hmm/{sample}_diagram.txt'
    params:
        extra = '',
    log:
        "logs/cnvkit/hmm/diagram/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'

rule cnvkit_heatmap_hmm:
    input:
        segments = expand("results/cnvkit/hmm/{sample}.cns", sample=samples.index)
    output:
        'results/cnvkit/hmm/heatmap.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/hmm/heatmap.txt'
    log:
        "logs/cnvkit/hmm/heatmap.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py heatmap -o {output} {input.segments} &> {log}'

rule export_seg_hmm:
    # Export the segmentation in DNAcopy format, i.e. create .seg file
    input:
        cns = 'results/cnvkit/hmm/{sample}.cns',
    output:
        'results/cnvkit/hmm/{sample}.seg'
    benchmark:
        'benchmarks/cnvkit/hmm/export_seg/{sample}.txt'
    params:
        extra = '--enumerate-chroms',
    log:
        "logs/cnvkit/hmm/export_seg/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py export seg {input.cns} -o {output} {params.extra} 2> {log}'

rule cnvkit_call_hmm:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output:
        'results/cnvkit/hmm/{sample}.call.cns'
    benchmark:
        'benchmarks/cnvkit/hmm/{sample}.call.txt'
    params:
        extra = '',
    log:
        "logs/cnvkit/hmm/call/{sample}.log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py call {input.copy_ratio} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'