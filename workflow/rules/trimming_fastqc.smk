rule fastqc_raw_reads:
  input:
    lambda wildcards: samples.at[wildcards.sample, 'fq1'] if wildcards.i == "1" else samples.at[wildcards.sample, 'fq2']
  output:
    html = "results/qc/fastqc/{sample}_{i}.html",
    zip = "results/qc/fastqc/{sample}_{i}_fastqc.zip"
  benchmark: "benchmarks/fastqc/{sample}_{i}.txt"
  log:
    "logs/fastqc/{sample}_{i}.log"
  conda:
    "../envs/qc.yaml"
  threads: 2 # two files per job
  wrapper:
    "v1.4.0/bio/fastqc"


rule trim:
  input:
    mate1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
    mate2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
  output:
    r1 = "results/trimmed/{sample}_1P.fq.gz",
    r2 = "results/trimmed/{sample}_2P.fq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired = "results/trimmed/{sample}_1U.fq.gz",
    r2_unpaired = "results/trimmed/{sample}_2U.fq.gz"
  benchmark: "benchmarks/trimmomatic/{sample}.txt"
  log:
    "logs/trimmomatic/{sample}.log"
  params:
    base = "results/trimmed/{sample}.fq.gz",
    # base = lambda wildcards: f"results/trimmed/{wildcards.sample}.fq.gz",
    adapter = config["adapter_sequences"]
  threads: 8
  resources:
    mem=lambda wildcards, attempt: '%dg' % (8 * attempt),
    slurm_partition = lambda wildcards, attempt: 'medium' if attempt > 1 else 'short',
    runtime=lambda wildcards, attempt: 24*60 if attempt > 1 else 4*60, # 4h=short partition limit, or 24h
    cores=lambda wc, threads: threads
  conda:
    "~/work/miniconda/envs/trim_map"
  shell:
    "trimmomatic PE "
      "-threads {threads} "
      "{input.mate1} {input.mate2} "
      "-baseout {params.base} "
      "ILLUMINACLIP:{params.adapter}:2:30:10 "
      "2> {log}"


rule fastqc_trimmed_reads:
  input:
    "results/trimmed/{sample}_{i}P.fq.gz"
  output:
    html = "results/qc/fastqc_trim/{sample}_{i}P.html",
    zip = "results/qc/fastqc_trim/{sample}_{i}P_fastqc.zip"
  benchmark: "benchmarks/fastqc_trim/{sample}_{i}P.txt"
  log:
    "logs/fastqc_trim/{sample}_{i}P.log"
  conda:
    "../envs/qc.yaml"
  threads: 2 # two files per job
  wrapper:
    "v1.4.0/bio/fastqc"

rule multiqc_fastqc_only:
  input:
    fastqc_html = expand("results/qc/fastqc/{sample}_{i}.html", sample = list(samples.index), i = ["1", "2"]),
    fastqc_trim_html = expand("results/qc/fastqc_trim/{sample}_{i}P.html", sample = list(samples.index), i = ["1", "2"]),
  output:
    "results/qc/multiqc_report_fastqc.html"
  benchmark: "benchmarks/multiqc_fastqc.txt"
  priority: -2 # error-prone rule, run it last
  log:
    "logs/multiqc.log"
  resources:
    mem=lambda wildcards, attempt: '%dG' % (16 * attempt),
    runtime=24*60, # 24h
    slurm_partition='medium'
  conda:
    "../envs/qc.yaml"
  shell:
    "multiqc ./results/qc/ -o results/qc 2> {log}" 
    "multiqc --force -o results/qc -n multiqc_report_fastqc.html results/qc/fastqc_trim results/qc/fastqc"

