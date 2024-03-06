rule fastq_to_bam:
    input:
        mate1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        mate2 = lambda wildcards: samples.at[wildcards.sample, 'fq2'],
    output:
        "results/umi_mapping/{sample}_1_fastqtobam.bam"
    conda:
        "../envs/umi_map.yaml"
    params:
        assay = "HS2-HRD"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        cores=lambda wc, threads: threads
    log: "logs/fastq_to_bam/{sample}.log",
    shell:
        """
        fgbio --tmp-dir=$TMPDIR \
            FastqToBam \
            -i {input.mate1} {input.mate2} \
            -o {output} \
            --sample={wildcards.sample} \
            --library={params.assay} \
            --sort=true \
            -Xmx{resources.mem} 2> {log}
        """


rule extract_umis_from_bam:
    input:
        "results/umi_mapping/{sample}_1_fastqtobam.bam",
    output:
        "results/umi_mapping/{sample}_2_extract_umis.bam"
    params:
        readstruc = "3M2S75T"
    conda:
        "../envs/umi_map.yaml"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        cores=lambda wc, threads: threads
    log: "logs/extract_umis/{sample}.log",
    shell:
        """
        fgbio --tmp-dir=$TMPDIR \
            ExtractUmisFromBam \
            -i {input} \
            -o {output} \
            -r {params.readstruc} {params.readstruc} \
            -t ZA ZB \
            -s RX \
            -Xmx{resources.mem} 2> {log}
        """


rule picard_samtofastq:
    input:
        "results/umi_mapping/{sample}_2_extract_umis.bam"
    output:
        "results/umi_mapping/{sample}_3_samtofastq.fastq"
    conda:
        "../envs/umi_map.yaml"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        cores=lambda wc, threads: threads
    log: "logs/samtofastq/{sample}.log",
    shell:
        "picard SamToFastq I={input} F={output} INTERLEAVE=TRUE -Xmx{resources.mem} 2> {log}"


rule bwa_mem_samples_umi:
    input:
        ref_index=config['reference'+'_'+config['genome_version']]["index"],
        reads="results/umi_mapping/{sample}_3_samtofastq.fastq",
        idx=multiext(stem, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output: "results/umi_mapping/{sample}_4_bwamem.bam"
    benchmark: 'benchmarks/bwa_mem_umi_mapping/{sample}.txt'
    log: "logs/bwa_mem/{sample}.log",
    threads: 16
    params:
        stem = stem,
        extra = "-p", # secondary alignment flags will cause MergeBamAlignment to fail
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        runtime=24*60, # 24h
        slurm_partition='medium',
        cores=lambda wc, threads: threads
    conda: "../envs/umi_map.yaml"
    log: "logs/bwa_mem/{sample}.log"
    shell: "bwa mem {params.extra} -t {threads} {params.stem} {input.reads} | samtools view -bS -@ {threads} 2> {log} > {output}"


rule picard_sortsam:
    input:
        "results/umi_mapping/{sample}_4_bwamem.bam"
    output:
        "results/umi_mapping/{sample}_5_sortsam.bam"
    resources:
        mem=lambda wildcards, attempt: '%dg' % (12 * attempt)
    conda:
        "../envs/umi_map.yaml"
    threads: 8
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        cores=lambda wc, threads: threads
    log: "logs/sortsam/{sample}.log"
    shell:
        """picard SortSam \
            -XX:ParallelGCThreads={threads} \
            I={input} \
            O={output} \
            SO=coordinate \
            TMP_DIR=$TMPDIR \
            -Xmx{resources.mem} 2> {log}
        """


rule MergeBamAlignment:
    input:
        unmapped = "results/umi_mapping/{sample}_2_extract_umis.bam",
        mapped = "results/umi_mapping/{sample}_5_sortsam.bam",
        genome=config['reference'+'_'+config['genome_version']]["fasta"]
    output:
        "results/umi_mapping/{sample}_6_mergebamalignment.bam"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (32 * attempt^2),
        slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
        runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60, # 4h=short partition limit, or 24h
        cores=lambda wc, threads: threads
    conda:
        "../envs/umi_map.yaml"
    log: "logs/MergeBamAlignment/{sample}.log"
    shell:
        "picard MergeBamAlignment "
            "-Xmx{resources.mem} "
            "UNMAPPED={input.unmapped} "
            "ALIGNED={input.mapped} "
            "O={output} "
            "R={input.genome} "
            "SO=coordinate "
            "ALIGNER_PROPER_PAIR_FLAGS=true "
            "MAX_GAPS=1 "
            "ORIENTATIONS=FR "
            "VALIDATION_STRINGENCY=SILENT "
            "CREATE_INDEX=true "
            "TMP_DIR=$TMPDIR "
            "MAX_RECORDS_IN_RAM=2000000  2> {log}"


rule UmiAwareMarkDuplicatesWithMateCigar:
    input:
        "results/umi_mapping/{sample}_6_mergebamalignment.bam"
    output:
        # FOLDER + "/" + PREFIX + ".mapped.annotated_umi.rmpcrduplex.bam",
        # FOLDER + "/" + PREFIX + "_ngspipmetrics.txt",
        # FOLDER + "/" + PREFIX + "_ngspipumimetrics.txt"
        rmpcrduplex = "results/umi_mapping/{sample}_7_UmiAwareMarkDuplicatesWithMateCigar.bam",
        ngspipmetrics = "results/umi_mapping/{sample}.ngspipmetrics.txt",
        ngspipumimetrics = "results/umi_mapping/{sample}.ngspipumimetrics.txt"
    conda:
        "../envs/umi_map.yaml"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt),
        cores=lambda wc, threads: threads
    log: "logs/UmiAwareMarkDuplicatesWithMateCigar/{sample}.log"
    shell:
        """
        picard UmiAwareMarkDuplicatesWithMateCigar \
            I={input} \
            O={output.rmpcrduplex} \
            M={output.ngspipmetrics} \
            UMI_METRICS={output.ngspipumimetrics} \
            DUPLEX_UMI=true \
            MAX_RECORDS_IN_RAM=2000000 \
            TMP_DIR=$TMPDIR \
            -Xmx{resources.mem} 2> {log}
        """


rule clipbam:
    input:
        rmpcrduplex = "results/umi_mapping/{sample}_7_UmiAwareMarkDuplicatesWithMateCigar.bam",
        genome=config['reference'+'_'+config['genome_version']]["fasta"]
    output:
        "results/umi_mapping/{sample}_8_clipbam.bam"
    conda:
        "../envs/umi_map.yaml"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: '%dg' % (16 * attempt^2),
        cores=lambda wc, threads: threads
    log: "logs/clipbam/{sample}.log"
    shell:
        """
        fgbio --tmp-dir=$TMPDIR \
            ClipBam \
            --input={input.rmpcrduplex} \
            --output={output} \
            --ref={input.genome} \
            --clipping-mode=Hard \
            --clip-overlapping-reads=true \
            -Xmx{resources.mem} 2> {log}
        """


rule index_clip_bam:
    input:
        "results/umi_mapping/{sample}_8_clipbam.bam"
    output:
        "results/umi_mapping/{sample}_8_clipbam.bam.bai"
    threads: 1
    resources:
        cores=lambda wc, threads: threads
    conda:
        "../envs/umi_map.yaml"
    log: "logs/index_clip_bam/{sample}.log"
    shell:
        "samtools index {input} 2> {log}"