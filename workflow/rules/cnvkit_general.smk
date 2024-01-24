import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# produced bed files
my_targets = 'results/cnvkit/general/my_targets.bed'
my_antitargets = 'results/cnvkit/general/my_antitargets.bed'

autobin_targets = 'results/cnvkit/general/'+bedname+'_target.bed'
autobin_antitargets = 'results/cnvkit/general/'+bedname+'_antitarget.bed'

bedname=bedname

# right bam files according to method
method_specific_bam = 'results/bam_sorted_bwa/{sample}_sorted_marked.bam' if not config["amplicon"] else 'results/bam_sorted_bwa/{sample}_sorted.bam'

# use remote file instead of cnvkit.py access output which causes problems
rule download_mappability:
    input:
        HTTP.remote("github.com/etal/cnvkit/raw/master/data/access-5k-mappable.hg19.bed", keep_local=True)
    output:
        config["mappability"]
    benchmark:
        "benchmarks/cnvkit/general/download_mappability.txt"
    log:
        "logs/cnvkit/general/download_mappability/log",
    shell:
        "mv {input} {output}"

rule cnvkit_target:
    input:
        config["bed_w_chr"],
    output:
        my_targets,
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

rule cnvkit_antitarget:
    input:
        bed = config["bed_w_chr"],
        access = config["mappability"],
    output:
        my_antitargets,
    benchmark:
        "benchmarks/cnvkit/general/antitarget.txt"
    log:
        "logs/cnvkit/general/antitarget/log",
    params:
        extra = '',
    conda:
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py antitarget {input.bed} --access {input.access} -o {output} {params.extra} 2> {log}'

rule cnvkit_autobin:
    input:
        bams = expand(method_specific_bam, sample=samples.index), # maybe use only fraction of samples here?
        targets = my_targets,
        access = config["mappability"],
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
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py autobin {input.bams} --targets {input.targets} --access {input.access} --target-output-bed {output.target} --antitarget-output-bed {output.antitarget} --method {params.method} {params.extra} 2> {log}'


rule cnvkit_coverage:
    input:
        bam = method_specific_bam,
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
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py coverage {input.bam} {input.targets} --processes {threads} -o {output.target_coverage} {params.extra} && '
        'cnvkit.py coverage {input.bam} {input.antitargets} --processes {threads} -o {output.antitarget_coverage} {params.extra} 2> {log}'

rule cnvkit_generic_ref:
    input:
        fasta=config["ref"],
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
        "../envs/primary_env.yaml"
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
        "../envs/primary_env.yaml"
    shell:
        'cnvkit.py fix {input.target_coverage} {input.antitarget_coverage} {input.reference} -o {output} {params.extra} {params.amplicon} 2> {log}'
