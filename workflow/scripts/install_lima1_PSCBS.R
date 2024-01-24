# Load the necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install the lima1/PSCBS package with the specified ref
BiocManager::install("remotes")
BiocManager::install("lima1/PSCBS", ref = "add_dnacopy_weighting")
