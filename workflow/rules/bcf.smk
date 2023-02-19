rule bcf_stats:
    input:
        #"{prefix}" # input: bvf/vcf file
        "results/bcf/call/{sample}.calls.bcf"
    output:
        #"{prefix}.stats.txt"
        "results/bcf/stats/{sample}.stats.txt"
    log:
        "logs/bcf_stats/{sample}.bcftools.stats.log"
    params:
        ""
    conda:
        "../envs/bcf.yaml"
    wrapper:
        "v1.7.1/bio/bcftools/stats"

rule bcftools_call:
    input:
        pileup="results/bcf/pileup/{sample}.pileup.bcf",
    output:
        calls="results/bcf/call/{sample}.calls.bcf",
    params:
        #caller="-m", # valid options include -c/--consensus-caller or -m/--multiallelic-caller
        caller = config["caller"],
        #options="--ploidy 1 --prior 0.001",
        options = config["caller_options"]
    log:
        "logs/bcftools_call/{sample}.log",
    conda:
        "../envs/bcf.yaml"
    wrapper:
        "v1.7.1/bio/bcftools/call"

rule bcftools_mpileup:
    input:
        index=config["ref_index"],
        ref=config["ref"], # this can be left out if --no-reference is in options
        alignments="results/mapped/{sample}.bam",
    output:
        pileup="results/bcf/pileup/{sample}.pileup.bcf",
    params:
        #options="--max-depth 100 --min-BQ 15",
        options = config["mpileup_options"],
    log:
        "logs/bcftools_mpileup/{sample}.log",
    conda:
        "../envs/bcf.yaml"
    threads: 
        8
    shell:
        "bcftools mpileup {params.options} --threads {threads} -o {output.pileup} -f {input.ref} {input.alignments} 2> {log}"
    #wrapper:
        #"v1.7.1/bio/bcftools/mpileup"

rule index:
    input:
        "results/bam_sorted_bwa/{sample}_sorted.bam"
    output:
        "results/bam_sorted_bwa/{sample}_sorted.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    threads:
        4
    conda:
        "../envs/stats.yaml"
    shell:
        "samtools index {input} -@ {threads} 2> {log}"
