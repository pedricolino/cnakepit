bedtargets = 'results/cnvkit/general/'+bedname+'_target.bed'
bedantitargets = 'results/cnvkit/general/'+bedname+'_antitarget.bed'
amplicontargets = 'results/cnvkit/general/'+bedname+'_amplicon.bed'
bedname=bedname

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# use remote file instead of cnvkit.py access output which causes problems
rule get_mappability:
    input:
        HTTP.remote("github.com/etal/cnvkit/raw/master/data/access-5k-mappable.hg19.bed", keep_local=True)
    output:
        config["mappability"]
    benchmark:
        "benchmarks/cnvkit/general/get_mappability.txt"
    log:
        "logs/cnvkit/general/get_mappability/log",
    shell:
        "mv {input} {output}"

# only for amplicon sequencing
if config["amplicon"]:
    rule cnvkit_target:
        input:
            config["bed_w_chr"],
        output:
            amplicontargets,
        benchmark:
            "benchmarks/cnvkit/general/target.txt"
        log:
            "logs/cnvkit/general/target/log",
        params:
            extra = '--split',
        conda:
            "../envs/primary_env.yaml"
        shell:
            'cnvkit.py target {input} -o {output} {params.extra} 2> {log}'

rule cnvkit_autobin:
    input:
        bams = expand("results/bam_sorted_bwa/{sample}_sorted.bam", sample=samples.index), # maybe use only fraction of samples here?
        targets = config["bed_w_chr"] if not config["amplicon"] else amplicontargets,
        access = config["mappability"],
    output:
        target = bedtargets,
        antitarget = bedantitargets,
    resources:
        mem_mb=8000 # otherwise 4TB for 600 samples with default Snakemake...
    benchmark:
        "benchmarks/cnvkit/general/autobin.txt"
    params:
        extra = '--method amplicon',
        samplenames = samples.index
    log:
        "logs/cnvkit/general/autobin/log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py autobin {input.bams} --targets {input.targets} --access {input.access} --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} {params.extra} 2> {log}'


rule cnvkit_coverage:
    input:
        bam = 'results/bam_sorted_bwa/{sample}_sorted.bam',
        targets = bedtargets if not config["amplicon"] else amplicontargets,
        antitargets = bedantitargets,
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
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py coverage {input.bam} {input.targets} --processes {threads} -o {output.target_coverage} {params.extra} && '
        'cnvkit.py coverage {input.bam} {input.antitargets} --processes {threads} -o {output.antitarget_coverage} {params.extra} 2> {log}'

rule cnvkit_ref_generic:
    input:
        fasta=config["ref"],
        targets = bedtargets if not config["amplicon"] else amplicontargets,
    output:
        FlatReference_cnn = 'results/cnvkit/general/FlatReference.cnn',
    benchmark: "benchmarks/cnvkit/general/ref_generic.txt"
    params:
        extra = '--sample-sex female',
        amplicon = '--no-edge' if config["amplicon"] else '',
    log:
        "logs/cnvkit/general/ref_generic/log",
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py reference -o {output.FlatReference_cnn} -f {input.fasta} -t {input.targets} {params.extra} {params.amplicon} 2> {log}'

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
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py fix {input.target_coverage} {input.antitarget_coverage} {input.reference} -o {output} {params.extra} {params.amplicon} 2> {log}'
