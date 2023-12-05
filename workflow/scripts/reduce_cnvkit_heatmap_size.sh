#!/bin/bash

# Downloaded script is licensed under the BSD 3-clause license, see script source code.
shrinkpdf_script="workflow/scripts/shrinkpdf.sh"

if [ ! -f "$shrinkpdf_script" ]; then
    curl -o "$shrinkpdf_script" https://raw.githubusercontent.com/aklomp/shrinkpdf/master/shrinkpdf.sh
    chmod +x "$shrinkpdf_script"
fi

# find heatmap pdfs in results/cnvkit subfolders
find results/cnvkit -name "heatmap.cnv.pdf" -print0 | while read -d $'\0' file
do
    outfile="${file//heatmap/heatmap_compressed}"
    echo "Reducing size of $file and save it to $outfile"
    "./$shrinkpdf_script" -o "$outfile" "$file"
done

