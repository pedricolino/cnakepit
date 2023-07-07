#ref_ref = str(Path("results") / "reference" / ref_file.stem)
stem = str(Path("results") / "ref" / ref_file.stem)


if config["DNA_seq"]:
    # bwa if DNA sequencing
    rule bwa_index:
        input:
            ref = ref_file,
            #ref = "results/filter_ref/hg38_filtered.fa",
        output:
            #idx=multiext(ref_ref, ".amb", ".ann", ".bwt", ".pac", ".sa"),
            idx=multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        log:
            "logs/bwa_index/bwa_index.log",
        params:
            algorithm="bwtsw",
        conda:
            "../envs/map.yaml"
        wrapper:
            "v1.7.0/bio/bwa/index"

    rule bwa_mem:
        input:
            reads=["results/trimmed/{sample}_1P.fq.gz", "results/trimmed/{sample}_2P.fq.gz"],
            idx=multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            "results/mapped/{sample}.bam",
        log:
            "logs/bwa_mem/{sample}.log",
        params:
            extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
            sorting="none",  # Can be 'none', 'samtools' or 'picard'.
            sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
            sort_extra="",  # Extra args for samtools/picard.
        threads: 16
        conda:
            "../envs/map.yaml"
        wrapper:
            "v1.7.0/bio/bwa/mem"

    ### sort and index mapped files so they can be imported to IGV
    #rule sort_bwa:
    #    input:
    #        "results/mapped/{sample}.bam"
    #    output:
    #        "results/bam_sorted_bwa/{sample}_sorted.bam"
    #    log:
    #        "logs/samtools/sort_bwa/{sample}.log"
    #    threads:
    #       8
    #    conda:
    #        "../envs/qc_map.yaml"
    #    shell:
    #        "samtools sort -o {output} {input} -@ {threads} 2> {log}"

    rule index_bwa:
        input:
            "results/bam_sorted_bwa/{sample}_sorted.bam"
        output:
            "results/bam_sorted_bwa/{sample}_sorted.bam.bai"
        log:
            "logs/samtools/index_bwa/{sample}.log"
        threads: 8
        conda:
            "../envs/qc_map.yaml"
        params:
            extra="",  # optional params string
        wrapper:
            "v2.0.0/bio/samtools/index"

else:
    # STAR if RNA sequencing
    rule star_pe_multi:
        input:
            # use a list for multiple fastq files for one sample
            # usually technical replicates across lanes/flowcells
            fq1=["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
            # paired end reads needs to be ordered so each item in the two lists match
            fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
            # path to STAR reference genome index
            idx="index",
        output:
            # see STAR manual for additional output files
            aln="star/pe/{sample}/pe_aligned.sam",
            log="logs/pe/{sample}/Log.out",
            sj="star/pe/{sample}/SJ.out.tab",
        log:
            "logs/star/{sample}.log",
        params:
            # optional parameters
            extra="",
        threads: 8
        conda:
            "../envs/map_rna.yaml" 
        wrapper:
            "v1.14.1/bio/star/align"

    rule star_index:
        input:
            fasta="{genome}.fasta",
        output:
            directory("{genome}"),
        threads: 1
        params:
            extra="",
        log:
            "logs/star/star_index/star_index_{genome}.log",
        conda:
            "../envs/map_rna.yaml" 
        wrapper:
            "v1.14.1/bio/star/index"
