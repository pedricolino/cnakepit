# Read the input TSV file
input_file <- if (length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "all_samples_filepath_size.tsv"
}

library(ggplot2)
df <- read.table(input_file, sep = "\t", header = TRUE)

# Save plot as PNG
png(
    filename = "file_size_distribution.png",
    width = 480, height = 770
)

cutoff_filesize <- 1e8

ggplot(df, aes(x = log10(FileSize))) +
    geom_histogram(bins = 100, fill = "skyblue", color = "black") +
    labs(
        x = "File size)",
        y = "Number of samples",
    ) +
    scale_x_continuous(
        breaks = 1:11,
        labels = c("10 B", "100 B", "1 KB", "10 KB", "100 KB", "1 MB", "10 MB", "100 MB", "1 GB", "10 GB", "100 GB"),
        limits = c(1, 11) # Add this line
    ) +
    geom_vline(xintercept = log10(cutoff_filesize), linetype = "dashed", color = "red")

dev.off()
