#!/bin/bash

Rscript workflow/scripts/collect_purecn_optima.R $(ls -d results/purecn/*/ | xargs)