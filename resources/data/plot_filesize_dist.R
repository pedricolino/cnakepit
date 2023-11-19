# Read the input TSV file
input_file <- if (length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "all_samples_filepath_size.tsv"
}

df <- read.table(input_file, sep="\t", header=TRUE)

# Save plot as PNG
png(filename="file_size_distribution.png",
    width = 480, height = 770)

# Plot file size in log10 scale
hist(log10(df$FileSize), xlab="File size", main="File Size Distribution", breaks=50, xaxt='n', col="lightblue")
# Set x-axis labels to file size units
axis(side = 1, at = 1:11, labels=c("10 B","100 B","1 KB","10 KB","100 KB","1 MB","10 MB","100 MB","1 GB","10 GB", "100 GB"))

# Save the plot
dev.off()

