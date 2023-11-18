#!/bin/bash

# Visualize the directed acyclic graph (DAG) of the workflow

dir=workflow/DAG/
mkdir -p $dir

# create dot file to customize output later
snakemake --rulegraph > ${dir}DAG.dot

# create a provisional png file
dot -Tpng ${dir}DAG.dot > ${dir}DAG.png