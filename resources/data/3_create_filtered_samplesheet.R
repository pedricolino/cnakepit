# Read the input TSV file
input_file <- if (length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "all_samples_filepath_size.tsv"
}

df <- read.table(input_file, sep="\t", header=TRUE)

# Filter out samples with file size less than 100000000 bytes
filtered_df <- df[df$FileSize >= 100000000, ]

# Group by sample name and get the first and second file paths
df_grouped <- split(filtered_df$FilePath, filtered_df$SampleName)
df_grouped <- lapply(df_grouped, function(x) c(x[1], if(length(x) > 1) x[2] else NA))
df_grouped <- data.frame(sample=names(df_grouped), do.call(rbind, df_grouped))
names(df_grouped) <- c("sample", "fq1", "fq2")

# Save the new DataFrame as a TSV file
write.table(df_grouped, file="sample_sheet.tsv", sep="\t", row.names=FALSE)

# Filtered out samples
filtered_samples <- df[!df$SampleName %in% df_grouped$sample, ]
write.table(filtered_samples, file="filtered_samples.tsv", sep="\t", row.names=FALSE)
# Filter out samples with file size less than 100000000 bytes
filtered_df <- df[df$FileSize >= 100000000, ]

# Group by sample name and get the first and second file paths
df_grouped <- split(filtered_df$FilePath, filtered_df$SampleName)
df_grouped <- lapply(df_grouped, function(x) c(x[1], if(length(x) > 1) x[2] else NA))
df_grouped <- data.frame(sample=names(df_grouped), do.call(rbind, df_grouped))
names(df_grouped) <- c("sample", "fq1", "fq2")

# Save the new DataFrame as a TSV file
write.table(df_grouped, file="sample_sheet.tsv", sep="\t", row.names=FALSE)

# Filtered out samples
filtered_samples <- df[!df$SampleName %in% df_grouped$sample, ]
write.table(filtered_samples, file="filtered_samples.tsv", sep="\t", row.names=FALSE)