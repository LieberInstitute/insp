#!/bin/bash
#$ -cwd
#$ -N NC1_rse
#$ -o ./NC1_rse.txt
#$ -e ./NC1_rse.txt
#$ -m e

## Usage:
# qsub NC1_rse.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.3.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/ajaffe/data/lab/insp/NC1/NC1.Rdata" -s "sumsNC1" -c 1
