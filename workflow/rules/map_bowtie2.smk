# optional mapper bowtie: 

rule bowtie2_build:
    input:
        ref = ref_file,
    output:
        multiext(
            "results/bowtie_index/hs38_genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra="",  # optional parameters
    threads: 8
    conda:
        env_prefix + 'map_bowtie2' + env_suffix
    wrapper:
        "v1.7.0/bio/bowtie2/build"

rule bowtie2:
    input:
        sample=["results/trimmed/{sample}_1P.fq.gz", "results/trimmed/{sample}_2P.fq.gz"],
        idx=multiext(
            "results/bowtie_index/hs38_genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "results/bowtie_mapped/{sample}.bam",
    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    conda:
        env_prefix + 'map_bowtie2' + env_suffix
    wrapper:
        "v1.7.0/bio/bowtie2/align"