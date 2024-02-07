import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# produced bed files
my_targets = 'results/cnvkit/general/my_targets.bed'
my_antitargets = 'results/cnvkit/general/my_antitargets.bed'

autobin_targets = 'results/cnvkit/general/'+bedname+'_target.bed'
autobin_antitargets = 'results/cnvkit/general/'+bedname+'_antitarget.bed'

bedname=bedname

# use remote file instead of cnvkit.py access output which causes problems
rule download_mappability:
    input:
        HTTP.remote(config["mappability"]["link"], keep_local=True)
    output:
        config["mappability"]["bed"]
    benchmark:
        "benchmarks/cnvkit/general/download_mappability.txt"
    log:
        "logs/cnvkit/general/download_mappability/log",
    shell:
        "mv {input} {output}"

rule download_sv_blacklist:
    input:
        HTTP.remote(config["sv_blacklist"]["link"], keep_local=True)
    output:
        config["sv_blacklist"]["bed"]
    benchmark:
        "benchmarks/cnvkit/general/download_sv_blacklist.txt"
    log:
        "logs/cnvkit/general/download_sv_blacklist/log",
    shell:
        "mv {input} {output}"

rule cnvkit_access:
    input:
        ref = config["reference"]["fasta"],
        sv_blacklist = config["sv_blacklist"]["bed"],
    output:
        'results/cnvkit/general/access.bed'
    benchmark:
        "benchmarks/cnvkit/general/access.txt"
    log:
        "logs/cnvkit/general/access/log",
    params:
        extra = '',
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py access {input.ref} --exclude {input.sv_blacklist} -o {output} {params.extra} 2> {log}'

rule cnvkit_target:
    input:
        config["panel_design"],
    output:
        my_targets,
    benchmark:
        "benchmarks/cnvkit/general/target.txt"
    log:
        "logs/cnvkit/general/target/log",
    params:
        extra = '--split',
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py target {input} -o {output} {params.extra} 2> {log}'

rule cnvkit_antitarget:
    input:
        bed = config["panel_design"],
        access = mappability
    output:
        my_antitargets,
    benchmark:
        "benchmarks/cnvkit/general/antitarget.txt"
    log:
        "logs/cnvkit/general/antitarget/log",
    params:
        extra = '',
    conda:
        "../envs/cnv_calling.yaml"
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
        "benchmarks/cnvkit/general/autobin.txt"
    params:
        extra = '',
        method = 'amplicon' if config["amplicon"] else 'hybrid',
        samplenames = samples.index
    log:
        "logs/cnvkit/general/autobin/log",
    conda:
        "../envs/cnv_calling.yaml"
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
    benchmark: "benchmarks/cnvkit/general/coverage/{sample}.txt"
    params:
        extra = '',
    threads: 8
    log:
        "logs/cnvkit/general/coverage/{sample}.log",
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py coverage {input.bam} {input.targets} --processes {threads} -o {output.target_coverage} {params.extra} && '
        'cnvkit.py coverage {input.bam} {input.antitargets} --processes {threads} -o {output.antitarget_coverage} {params.extra} 2> {log}'

rule cnvkit_generic_ref:
    input:
        fasta=config["reference"]["fasta"],
        targets = autobin_targets,
        antitargets = autobin_antitargets,
    output:
        FlatReference_cnn = 'results/cnvkit/general/FlatReference.cnn',
    benchmark: "benchmarks/cnvkit/general/ref_generic.txt"
    params:
        extra = '--sample-sex female',
        amplicon = '--no-edge' if config["amplicon"] else '',
    log:
        "logs/cnvkit/general/ref_generic/log",
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py reference -o {output.FlatReference_cnn} -f {input.fasta} -t {input.targets} -a {input.antitargets} {params.extra} {params.amplicon} 2> {log}'

rule cnvkit_fix:
    input:
        target_coverage = 'results/cnvkit/general/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/general/{sample}.antitargetcoverage.cnn',
        reference = 'results/cnvkit/general/FlatReference.cnn',
    output:
        'results/cnvkit/general/{sample}.cnr'
    benchmark: "benchmarks/cnvkit/general/fix/{sample}.txt"
    params:
        extra = '',
        amplicon = '--no-edge' if config["amplicon"] else '',
    log:
        "logs/cnvkit/general/fix/{sample}.log",
    conda:
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py fix {input.target_coverage} {input.antitarget_coverage} {input.reference} -o {output} {params.extra} {params.amplicon} 2> {log}'


# CNVkit segmentation and downstream rules

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
    resources:
        mem=lambda wildcards, attempt: '%dG' % (4 * attempt),
    log:
        "logs/cnvkit/cbs/segment/{sample}.log",
    conda:
        "../envs/cnv_calling.yaml"
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
        "../envs/cnv_calling.yaml"
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
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py diagram {input.copy_ratio} --segment {input.segment} -o {output} {params.extra} 2> {log}'


rule cnvkit_heatmap_cbs:
    input:
        segments = expand("results/cnvkit/cbs/{sample}.cns", sample=samples.index)
    output:
        'results/cnvkit/cbs/heatmap.cnv.pdf'
    benchmark:
        'benchmarks/cnvkit/cbs/heatmap.txt'
    resources:
        mem=lambda wildcards, attempt: '%dG' % (32 * attempt),
        slurm_partition = 'medium'
    log:
        "logs/cnvkit/cbs/heatmap.log",
    conda:
        "../envs/cnv_calling.yaml"
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
        "../envs/cnv_calling.yaml"
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
        "../envs/cnv_calling.yaml"
    shell:
        'cnvkit.py call {input.cns} --vcf {input.germline_vcf} -o {output} {params.extra} 2> {log}'


# Do the same with the HMM segmentation algorithm

use rule cnvkit_segment_cbs as cnvkit_segment_hmm with:
    output: 'results/cnvkit/hmm/{sample}.cns',
    benchmark: 'benchmarks/cnvkit/hmm_segment/{sample}.txt'
    params: extra = '-m hmm --drop-low-coverage',
    log: "logs/cnvkit/hmm/segment/{sample}.log",

use rule cnvkit_scatter_cbs as cnvkit_scatter_hmm with:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output: 'results/cnvkit/hmm/{sample}_scatter.cnv.pdf'
    benchmark: 'benchmarks/cnvkit/hmm/{sample}_scatter.txt'
    log: "logs/cnvkit/hmm/scatter/{sample}.log",

use rule cnvkit_diagram_cbs as cnvkit_diagram_hmm with:
    input:
        copy_ratio = 'results/cnvkit/general/{sample}.cnr',
        segment = 'results/cnvkit/hmm/{sample}.cns',
    output: 'results/cnvkit/hmm/{sample}_diagram.cnv.pdf'
    benchmark: 'benchmarks/cnvkit/hmm/{sample}_diagram.txt'
    log: "logs/cnvkit/hmm/diagram/{sample}.log",

use rule cnvkit_heatmap_cbs as cnvkit_heatmap_hmm with:
    input: segments = expand("results/cnvkit/hmm/{sample}.cns", sample=samples.index)
    output: 'results/cnvkit/hmm/heatmap.cnv.pdf'
    benchmark: 'benchmarks/cnvkit/hmm/heatmap.txt'

use rule export_seg_cbs as export_seg_hmm with:
    input: cns = 'results/cnvkit/hmm/{sample}.cns',
    output: 'results/cnvkit/hmm/{sample}.seg'
    benchmark: 'benchmarks/cnvkit/hmm/export_seg/{sample}.txt'
    log: "logs/cnvkit/hmm/export_seg/{sample}.log",

use rule cnvkit_call_cbs as cnvkit_call_hmm with:
    input:
        cns = 'results/cnvkit/hmm/{sample}.cns',
        germline_vcf = "results/mutect2/germline/{sample}_germline.vcf"
    output: 'results/cnvkit/hmm/{sample}.icns' # integer copy number
    benchmark: 'benchmarks/cnvkit/hmm/{sample}.call.txt'
    log: "logs/cnvkit/hmm/call/{sample}.log",