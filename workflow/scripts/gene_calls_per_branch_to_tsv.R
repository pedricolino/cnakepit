library(tidyverse)

dir <- "results/collection_of_specific_samples/results"
all_genes_file <- "resources/paneldesign/all_genes.csv"
list.dirs(path = dir, full.names = TRUE, recursive = FALSE) -> pons

# for every folder, read in all _genes.csv files, extract the gene and type columns, and merge into one table with gene as row names and C in the columns names according to the sample name contained in the file name. Then, for every branch, write the table to a file called branch+'_gene_calls.tsv'.

all_genes <- read_csv(all_genes_file, col_names = FALSE) # %>% pull(X1)

empty_row_tables <- c()
list <- list()
for (pon in pons) {
    list.dirs(path = pon, full.names = TRUE, recursive = FALSE) -> branches
    for (branch in branches) {
        list.files(path = branch, pattern = "_genes.csv", full.names = TRUE) -> files
        if (length(files) == 0) {
            next
        }
        list <- list()
        for (file in files) {
            print(file)
            read_csv(file, show_col_types = FALSE) %>%
                select(gene.symbol, type) %>%
                filter(
                    !grepl("chr[x0-9]*:.*", gene.symbol),
                    !grepl("Antitarget.*", gene.symbol),
                    !grepl("x:.*", gene.symbol)
                ) %>%
                right_join(all_genes, by = c("gene.symbol" = "X1")) -> foo
            str_replace(file, "_genes.csv", "") %>%
                str_replace(".*/", "") -> name
            # foo %>% select(type) -> foo
            colnames(foo) <- c("gene", name)

            if (foo %>% nrow() == 0) {
                empty_row_tables <- c(empty_row_tables, name)
            } else {
                list[[name]] <- foo
            }
        }
        # merge list to table
        res <- reduce(list, full_join, by = "gene")
        branch %>%
            gsub(".*/results/", "", .) %>%
            gsub("/", "_", .) -> branchname
        write_tsv(res, paste0(dir, "/", branchname, "_gene_calls.tsv"))
    }
}
