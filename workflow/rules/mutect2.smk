#ref_ref = str(Path("results") / "reference" / ref_file.stem)

rule mutect2_bam:
    input:
        fasta=config["ref"],
        map="results/bam_sorted_bwa/{sample}_sorted.bam",
        dict="resources/reference/hs38.dict",
        idx="results/bam_sorted_bwa/{sample}_sorted.bam.bai",
        targets=config["bed"]
    output:
        vcf="results/mutect2/variant_bam/{sample}.vcf.gz",
        #bam="results/mutect2/variant_bam/{sample}.bam",
    #message:
        #"Testing Mutect2 with {wildcards.sample}"
    params:
        extra="--germline-resource $GERMLINE_RESOURCE --genotype-germline-sites true --genotype-pon-sites true --interval-padding 75"
    threads: 8
    log:
        "logs/mutect2/{sample}.log",
    conda:
        "../envs/mutect2.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I -L {input.targets} {input.map} -O {output.vcf} {params.extra} 2> {log}"

rule mutect2_ref_dict:
    input:
        config["ref"]
    output:
        "resources/reference/hs38.dict"
    conda:
        "../envs/stats.yaml"
    log:
        "logs/mutect2_dict/mutect2_dict.log"
    shell:
        "samtools dict {input} > {output} 2> {log}"

rule filter_mutect_calls:
    input:
        vcf="results/mutect2/variant_bam/{sample}.vcf.gz",
        reference=config["ref"],
    output:
        vcf_filt="results/mutect2/filtered/{sample}_filtered.vcf.gz",
    shell:
        "gatk FilterMutectCalls -R {input.reference} -V {input.vcf} -O {output.vcf_filt} 2> logs/mutect2_filter/{wildcards.sample}.log"