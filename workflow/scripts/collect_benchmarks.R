suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

benchmark_dir <- 'benchmarks'
target_file <- 'results/stats/benchmarks.tsv'
dir.create('results/stats', recursive = TRUE, showWarnings = FALSE)

# Remove the target file if it exists
if (file.exists(target_file)) {
    file.remove(target_file)
}

# List to store all dataframes
dfs <- list()

# Traverse through all subdirectories and files
files <- list.files(benchmark_dir, recursive = TRUE, full.names = TRUE)

for (file_path in files) {
    # Skip the target file
    if (file_path == target_file) {
        next
    }
    # Read the file into a DataFrame
    df <- read_tsv(file_path, show_col_types = FALSE)
    # Add a new column with the file path, and remove the benchmark directory from the path
    df <- df %>% 
        mutate(file_path = gsub("benchmarks/", "", file_path),
               # extract everything in front the first slash or dot
               rule = gsub("([^/\\.]+).*", "\\1", file_path))
    # Append the DataFrame to the list
    dfs <- append(dfs, list(df))
}

# Concatenate all dataframes
df_final <- bind_rows(dfs)

# Write the final DataFrame to a TSV file
write_tsv(df_final, target_file)