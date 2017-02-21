#!/bin/bash
#$ -cwd
#$ -N IGH_rse
#$ -o ./IGH_rse.txt
#$ -e ./IGH_rse.txt
#$ -m e

## Usage:
# qsub IGH_rse.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.3.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/lieber/ajaffe/lab/insp/IGH/IGH.Rdata" -s "sumsIGH" -c 2
