rule fastqc:
  input:
    lambda wildcards: samples.at[wildcards.sample, 'fq1'] if wildcards.i == "1" else samples.at[wildcards.sample, 'fq2']
  output:
    html = "results/qc/fastqc/{sample}_{i}.html",
    zip = "results/qc/fastqc/{sample}_{i}_fastqc.zip"
  log:
    "logs/fastqc/{sample}_{i}.log"
  conda:
    "../envs/primary_env.yaml"
  threads: 16
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
  log:
    "logs/trimmomatic/{sample}.log"
  params:
    base = lambda wildcards: f"results/trimmed/{wildcards.sample}.fq.gz",
    adapter = config["adapter"]
  threads: 16
  conda:
    "../envs/primary_env.yaml"
  shell:
    "trimmomatic PE -threads {threads} {input.mate1} {input.mate2} -baseout {params.base} ILLUMINACLIP:{params.adapter}:2:30:10 2> {log}"


rule fastqc_trim:
  input:
    "results/trimmed/{sample}_{i}P.fq.gz"
  output:
    html = "results/qc/fastqc_trim/{sample}_{i}P.html",
    zip = "results/qc/fastqc_trim/{sample}_{i}P_fastqc.zip"
  log:
    "logs/fastqc_trim/{sample}_{i}P.log"
  conda:
    "../envs/primary_env.yaml"
  threads: 16
  wrapper:
    "v1.4.0/bio/fastqc"


rule multiqc_fastqc:
  input:
    fastqc_html = expand("results/qc/fastqc/{sample}_{i}.html", sample = list(samples.index), i = ["1", "2"]),
    fastqc_trim_html = expand("results/qc/fastqc_trim/{sample}_{i}P.html", sample = list(samples.index), i = ["1", "2"]),
  output:
    "results/qc/multiqc_report_fastqc.html"
  log:
    "logs/multiqc.log"
  conda:
    "../envs/primary_env.yaml"
  shell:
    "multiqc ./results/qc/ -o results/qc 2> {log}" 
    "multiqc --force -o results/qc -n multiqc_report_fastqc.html results/qc/fastqc_trim results/qc/fastqc"

