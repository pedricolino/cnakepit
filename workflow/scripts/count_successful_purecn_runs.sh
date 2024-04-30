#!/bin/bash

find results/ -name '*rds*' | grep -o "results/[a-zA-Z_]*/[a-zA-Z_]*/*[^0-9]*/" | uniq -c >results/stats/successful_purecn_runs.txt
#find results/ -name '*rds*' | grep -o "results/[a-zA-Z_]*/[a-zA-Z_]*" | uniq -c >results/stats/successful_purecn_runs.txt
