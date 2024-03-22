#!/bin/bash

find results/ -name '*rds*' | sed 's/\AS-.*//g' | uniq -c > results/stats/successful_purecn_runs.txt