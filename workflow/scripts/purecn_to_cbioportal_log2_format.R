args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript purecn_to_cbioportal_log2_format.R <input_file_1> <input_file_2> ... <input_file_n> <output_file>")
}

suppressPackageStartupMessages({
	library(tidyverse)
	library(org.Hs.eg.db) # for entrez id
})

read_csv(args[-1], show_col_types = F) %>%
	# filter out antitargets, empty genes and duplicate genes
	filter(!gene.symbol %>% str_detect('__')) %>%
	dplyr::select(Sampleid, gene.symbol, seg.mean) %>%
	# 1 column per sample
	pivot_wider(names_from = Sampleid, values_from = seg.mean) %>%
	# add entrez id
	mutate(Entrez_Gene_Id = mapIds(
		org.Hs.eg.db, 
		keys = gene.symbol, 
		column = "ENTREZID", 
		keytype = "ALIAS"), # look up synonyms if no entrez id, e.g. C1orf170 (no ID) = PERM1 (has ID)
		) %>%
	# filter out rows without entrez IDs, a lot of them are called LOCxxxxxx
	filter(!is.na(Entrez_Gene_Id)) %>%
	relocate(Entrez_Gene_Id, .after = gene.symbol) %>%
	dplyr::rename(Hugo_Symbol = gene.symbol) %>%
	write_tsv(., args[1])