#ref_ref = str(Path("results") / "reference" / ref_file.stem)

rule mutect2_bam:
    input:
        fasta=config["ref"],
        map="results/bam_sorted_bwa/{sample}_sorted.bam",
        dict="resources/reference/hs38.dict",
        idx="results/bam_sorted_bwa/{sample}_sorted.bam.bai"
    output:
        vcf="results/mutect2/variant_bam/{sample}.vcf",
        #bam="results/mutect2/variant_bam/{sample}.bam",
    #message:
        #"Testing Mutect2 with {wildcards.sample}"
    threads: 8
    log:
        "logs/mutect2/{sample}.log",
    conda:
        "../envs/mutect2.yaml"
    shell:
        "gatk Mutect2 -R {input.fasta} -I {input.map} -O {output.vcf} 2> {log}"

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

