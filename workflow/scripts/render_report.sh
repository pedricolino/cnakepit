#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate R

Rscript -e "rmarkdown::render('summarise_benchmarks.Rmd')"
