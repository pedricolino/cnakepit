args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
    stop("At least one argument must be supplied (directory with PureCN results in sample subfolders).\n", call. = FALSE)
}

for (dir in args) {
    if (!dir.exists(dir)) {
        stop("Directory ", dir, " does not exist.\n", call. = FALSE)
    }

    # Find all .rds files in the directory and its subdirectories
    rds_files <- list.files(dir, pattern = "\\.rds$", recursive = TRUE)

    # Initialize counters
    processed_files <- 0
    total_files <- length(rds_files)

    # extract values from each model from each file
    list <- list()
    for (sample in rds_files) {
        rds <- readRDS(paste(dir, sample, sep = ""))
        sample_name <- gsub(".rds", "", basename(sample))
        nb <- length(rds$results)
        df <- data.frame("sample" = rep(sample_name, nb))
        for (i in 1:nb) {
            df[i, "candidate.id"] <- rds$results[[i]]$candidate.id
            df[i, "purity"] <- rds$results[[i]]$purity
            df[i, "ploidy"] <- rds$results[[i]]$ploidy
            df[i, "total.ploidy"] <- rds$results[[i]]$total.ploidy
            df[i, "log.likelihood"] <- rds$results[[i]]$log.likelihood
            df[i, "total.log.likelihood"] <- rds$results[[i]]$total.log.likelihood
            df[i, "log.ratio.sdev"] <- rds$results[[i]]$log.ratio.sdev
            df[i, "fraction.subclonal"] <- rds$results[[i]]$fraction.subclonal
        }
        list[[sample_name]] <- df

        # Update progress
        processed_files <- processed_files + 1
        cat(paste("\rProcessed", processed_files, "out of", total_files, "files in directory", dir))
    }

    df <- do.call(rbind, list) # list to df

    file <- paste("results/stats/", gsub(".*/", "", gsub("/$", "", dir)), "_ploidy_purity_local_optima.tsv", sep = "")

    write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

    cat(paste("\nOutput file saved as:", file, "\n"))
}
