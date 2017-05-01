#!/bin/bash
#$ -cwd
#$ -N IGH_rse_windows
#$ -o ./IGH_rse_windows.txt
#$ -e ./IGH_rse_windows.txt
#$ -m e

## Usage:
# qsub IGH_rse_windows.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.3.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/lieber/ajaffe/lab/insp/IGHwindows/IGHwindows.Rdata" -s "sumsIGHwindows" -c 1
