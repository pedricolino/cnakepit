#!/bin/bash

output="mapped_to_roi_secondary_alignments.tsv"
touch $output

echo -e "Sample\tMapped_to_ROI(%)\tSecondary_alignments" >$output
for file in $(ls -R */*.html); do
    name=$(echo $file | sed -e "s/\/qualimap.*//g")
    map_roi=$(cat $file | sed -n '181p' | sed -e "s/.*  \///g" | sed -e "s/%<\/td>//g" | sed -e "s/\s//g")
    sec=$(cat $file | sed -n '154p' | sed -e "s/<td class=column2>//g" | sed -e "s/<\/td>//g" | sed -e "s/,//g")
    echo -e "$name\t$map_roi\t$sec" >>$output
done

echo "File written to $output"
