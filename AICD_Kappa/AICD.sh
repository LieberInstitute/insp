#!/bin/bash
#$ -cwd
#$ -N AICD
#$ -o ./AICD.txt
#$ -e ./AICD.txt
#$ -m e

## Usage:
# qsub AICD.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.3.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/ajaffe/data/lab/insp/AICD_Kappa/AICD.Rdata" -s "sumsAICD" -c 1
