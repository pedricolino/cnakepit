rule add_padding:
  input:
    bed = config['panel_design'],
    refgenome = config['reference'+'_'+config['genome_version']]['fa_genome']
  output: "resources/paneldesign/padded.bed",
  log: "logs/mask_ref/add_padding.log"
  params: padding = config['padding']
  conda: "../envs/bedtools.yaml"
  shell: "bedtools slop -i {input.bed} -g {input.refgenome} -b {params.padding} > {output} 2> {log}"

rule sort_regions:
  input: "resources/paneldesign/padded.bed"
  output: "resources/paneldesign/padded_sorted.bed"
  log: "logs/mask_ref/sort_regions.log"
  conda: "../envs/bedtools.yaml"
  shell: "bedtools sort -i {input} > {output} 2> {log}"

rule complement_regions:
  input:
    bed = "resources/paneldesign/padded_sorted.bed",
    ref_index = config['reference'+'_'+config['genome_version']]['index']
  output: "resources/paneldesign/complement.bed"
  log: "logs/mask_ref/complement_regions.log"
  conda: "../envs/bedtools.yaml"
  shell: "bedtools complement -i {input.bed} -g {input.ref_index} > {output} 2> {log}"

rule mask_reference:
  input:
    ref = config['reference'+'_'+config['genome_version']]['fasta'],
    mask = "resources/paneldesign/complement.bed"
  output: config['masked_ref']
  log: "logs/mask_ref/mask_reference.log"
  conda: "../envs/bedtools.yaml"
  shell: "bedtools maskfasta -fi {input.ref} -bed {input.mask} -fo {output} 2> {log}"
