#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
    stop("At least one argument must be supplied (directory with PureCN results in sample subfolders).\n", call. = FALSE)
}

for (dir in args) {
    if (!dir.exists(dir)) {
        stop("Directory ", dir, " does not exist.\n", call. = FALSE)
    }

    # for every folder, find the files ending with _genes.csv
    files <- list.files(dir, pattern = "\\_genes.csv$", recursive = TRUE)

    # Initialize counters
    processed_files <- 0
    total_files <- length(files)

    # initialize empty lists
    copy_numbers <- list()
    log2r <- list()

    for (sample in files) {
        # read in file
        table <- read.csv(paste(dir, "/", sample, sep = ""))

        # discard rows with Antitargets
        table_filtered <- table[which(!grepl("Antitarget.*", table$gene.symbol)), ]

        # discard rows without gene names
        table_filtered <- table_filtered[which(!grepl("chr[0-9]*:.*", table_filtered$gene.symbol)), ]

        # save the required columns
        absolute_copy_numbers <- table_filtered[, c("gene.symbol", "C")]
        log2_ratios <- table_filtered[, c("gene.symbol", "gene.mean")]

        # extract sample name and rename columns
        sample_name <- regmatches(sample, regexpr("AS-[0-9LR-]*", sample))
        colnames(absolute_copy_numbers) <- c("gene", sample_name)
        colnames(log2_ratios) <- c("gene", sample_name)

        # save tables in lists
        copy_numbers[[sample_name]] <- absolute_copy_numbers
        log2r[[sample_name]] <- log2_ratios

        # Update progress
        processed_files <- processed_files + 1
        cat(paste("\rProcessed", processed_files, "out of", total_files, "files in directory", dir))
    }

    cat("\nAll files haven been processed. Tables are now being prepared to be saved.")
    # lists to dataframes
    cns <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), copy_numbers)
    l2r <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), log2r)

    # choose output file names
    file_acn <- paste("results/stats/", gsub(".*/", "", gsub("/$", "", dir)), "_absolute_copy_numbers_all_samples.tsv", sep = "")
    file_l2r <- paste("results/stats/", gsub(".*/", "", gsub("/$", "", dir)), "_log2_ratios_all_samples.tsv", sep = "")

    # save tables in files
    cat("\nWriting tables as TSV files.")
    write.table(cns, file_acn, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
    write.table(l2r, file_l2r, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

    cat(paste("\nOutput files saved as", file_acn, "and", file_l2r, "\n"))
}
