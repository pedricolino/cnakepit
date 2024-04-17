#!/bin/bash

# Display Help
Help() {
   echo "This script finds filepaths and file size for filenames listed in the input file in a specified directory."
   echo
   echo -e "\033[34mSyntax: get_filepaths.sh [-h] inputFile [directory]\033[0m"
   echo "Options:  -h            Print this Help."
   echo "          directory     The directory to search in. Defaults to SIGN-OC directory if not provided."
   echo
}

# Process the input options. Add options as needed.
# If no options, run the main script.
if [ -z "$1" ]; then
   Help
   exit
else
   while getopts ":h" option; do
      case $option in
      h) # display Help
         Help
         exit
         ;;
      esac
   done
fi

# Directory to search in
dir=${1:-"/fast/work/groups/cubi/projects/2021-12-10_Keilholz_SIGN_OC/"}

all_files_file=all_files.txt

# Search for all fastq files in the directory
find $dir -name "*.fastq.gz" | grep "T1-DNA" >$all_files_file

# Output file
output_file="all_samples_filepath_size.tsv"

# Clear the output file and write the headers
echo -e "FileName\tSampleName\tFileSize\tFilePath" >$output_file

# Get total number of lines in the input file
total_lines=$(wc -l <"$all_files_file")

# Initialize a counter
counter=0

# Read the input file line by line
while IFS= read -r filepath; do
   # Increment the counter
   ((counter++))

   # Find the file in the directory
   filename=$(basename $filepath)

   # Extract the sample name
   sample_name=${filename%_R[12]*.fastq.gz}
   sample_name=$(sed -e "s/_Lx//g" <<<$sample_name)

   # Get the size of the file in bytes
   file_size=$(du -b $filepath | cut -f1)

   # Write the details to the output file
   echo -e "$filename\t$sample_name\t$file_size\t$filepath" >>$output_file

   # Display progress
   echo -ne "Processed $counter out of $total_lines rows.\r"

done <"$all_files_file"
