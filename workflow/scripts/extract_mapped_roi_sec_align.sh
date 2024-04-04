#!/bin/bash

output="results/stats/mapped_to_roi_secondary_alignments.tsv"
touch $output

echo -e "Sample\tMapped_to_ROI(%)\tSecondary_alignments" >$output
for file in $(ls -R results/qc/qualimap/qualimap_bwa/*/*.html); do
    name=$(echo $file | sed -e "s/.*qualimap_bwa.//g" | sed -e "s/.qualimapReport.html//g")
    map_roi=$(cat $file | sed -n '181p' | sed -e "s/.*  \///g" | sed -e "s/%<\/td>//g" | sed -e "s/\s//g")
    sec=$(cat $file | sed -n '154p' | sed -e "s/<td class=column2>//g" | sed -e "s/<\/td>//g" | sed -e "s/,//g")
    echo -e "$name\t$map_roi\t$sec" >>$output
done

echo "File written to $output"
