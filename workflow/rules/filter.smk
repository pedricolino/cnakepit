rule filter_ref:
  input:
    ref=config["ref"],
    bed=config["bed"],
  output:
    "results/filter_ref/hg38_filtered.fa"
  log:
    "logs/filter/filter_ref.log"
  threads: 8
  conda:
    "../envs/filter.yaml"
  shell:
    "bedtools getfasta -fi {input.ref} -bed {input.bed} -fo {output} 2> {log}"

rule idx_ref_f:
  input:
    "results/filter_ref/hg38_filtered.fa"
  output:
    "results/filter_ref/hg38_filtered.fa.fai"
  log:
    "logs/filter/idx_ref.log"
  threads: 8
  conda:
    "../envs/filter.yaml"
  shell:
    "samtools faidx {input} > {output} 2> {log}"


