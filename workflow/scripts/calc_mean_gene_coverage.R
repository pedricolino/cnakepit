folder_path <- "results/cnvkit/general"
file_list <- list.files(folder_path, pattern = "\\.targetcoverage\\.cnn$", full.names = TRUE)

print(paste("Number of detected files:", length(file_list)))

# calculate gene coverages, weighted by gene length
calc_mean_coverages <- function(file) {
    # read in the file
    df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    # calculate gene lengths from columns start and end
    df$length <- df$end - df$start + 1
    # calculate gene coverages from columns gene, start and end
    df$mean_coverage <- df$depth * df$length
    # remove comma and everything after it from column gene
    df$gene <- gsub(",.*", "", df$gene)
    # sum up the gene coverages and length per gene
    agg <- aggregate(cbind(mean_coverage, length) ~ gene, data = df, sum)
    # divide by gene length
    agg$mean_coverage <- agg$mean_coverage / agg$length
    # set gene column as row names
    rownames(agg) <- agg$gene
    # return the data frame
    agg[, c("length", "mean_coverage")]
}

# apply the function to all files
df_list <- lapply(file_list, calc_mean_coverages)

print("Mean coverages calculated.")

# combine the mean_coverage columns from all data frames into one data frame
df <- do.call(cbind, lapply(df_list, "[", "mean_coverage"))
# set the column names to the file names
colnames(df) <- gsub(".targetcoverage.cnn", "", gsub(".*/", "", file_list))
# row names to column
df$gene <- rownames(df)
# reorder the columns
df <- df[, c("gene", colnames(df)[-ncol(df)])]

# write the data frame to a file
out_folder <- "results/stats"
output <- paste(out_folder, "/gene_coverages_all_samples.tsv", sep = "")
write.table(df, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

print(paste("Mean coverages written to", output))
