#!/bin/bash
#$ -cwd
#$ -N IGH_windows
#$ -o ./IGH_windows.txt
#$ -e ./IGH_windows.txt
#$ -m e

## Usage:
# qsub IGH_windows.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.3.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/ajaffe/data/lab/insp/IGHwindows/IGHwindows.Rdata" -s "sumsIGHwindows" -c 1
