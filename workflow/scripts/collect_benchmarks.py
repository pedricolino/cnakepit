import os
import pandas as pd

benchmark_dir = 'benchmarks'
target_file = benchmark_dir+'/all_benchmarks.tsv'

# Remove the target file if it exists
if os.path.exists(target_file):
    os.remove(target_file)

# List to store all dataframes
dfs = []

# Traverse through all subdirectories and files
for root, dirs, files in os.walk(benchmark_dir):
    for file in files:
        # Construct full file path
        file_path = os.path.join(root, file)
        # Skip the target file
        if file_path == target_file:
            continue
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        # Add a new column with the file path, and remove the benchmark directory from the path
        df['file_path'] = file_path.replace("benchmarks/", "")
        # Append the DataFrame to the list
        dfs.append(df)

# Concatenate all dataframes
df_final = pd.concat(dfs, ignore_index=True)

# Write the final DataFrame to a TSV file
df_final.to_csv(target_file, sep='\t', index=False)