rule bcf_stats:
    input:
        #"{prefix}" # input: bvf/vcf file
        "results/bcf/call/{sample}.calls.bcf"
    output:
        #"{prefix}.stats.txt"
        "results/bcf/stats/{sample}.stats.txt"
    benchmark:
        "benchmarks/bcf_stats/{sample}.bcftools.stats.txt"
    log:
        "logs/bcf_stats/{sample}.bcftools.stats.log"
    params:
        ""
    conda:
        env_prefix + 'qc' + env_suffix
    wrapper:
        "v1.7.1/bio/bcftools/stats"

rule bcftools_call:
    input:
        pileup="results/bcf/pileup/{sample}.pileup.bcf",
    output:
        calls="results/bcf/call/{sample}.calls.bcf",
    benchmark:
        "benchmarks/bcf_call/{sample}.bcftools.call.txt"
    params:
        #caller="-m", # valid options include -c/--consensus-caller or -m/--multiallelic-caller
        caller = config["caller"],
        #options="--ploidy 1 --prior 0.001",
        options = config["caller_options"]
    log:
        "logs/bcftools_call/{sample}.log",
    conda:
        env_prefix + 'qc' + env_suffix
    wrapper:
        "v1.7.1/bio/bcftools/call"

rule bcftools_mpileup:
    input:
        index=config['reference'+'_'+config['genome_version']]["index"],
        ref=config['reference'+'_'+config['genome_version']]["fasta"], # this can be left out if --no-reference is in options
        alignments="results/mapped/{sample}.bam",
    output:
        pileup="results/bcf/pileup/{sample}.pileup.bcf",
    benchmark:
        "benchmarks/bcf_mpileup/{sample}.bcftools.mpileup.txt"
    params:
        #options="--max-depth 100 --min-BQ 15",
        options = config["mpileup_options"],
    log:
        "logs/bcftools_mpileup/{sample}.log",
    conda:
        env_prefix + 'qc' + env_suffix
    threads: 
        8
    shell:
        "bcftools mpileup {params.options} --threads {threads} -o {output.pileup} -f {input.ref} {input.alignments} 2> {log}"
    #wrapper:
        #"v1.7.1/bio/bcftools/mpileup"
