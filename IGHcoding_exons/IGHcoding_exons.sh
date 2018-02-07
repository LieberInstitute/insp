#!/bin/bash
#$ -cwd
#$ -N IGHcoding_exons
#$ -o ./IGHcoding_exons.txt
#$ -e ./IGHcoding_exons.txt
#$ -m e

## Usage:
# qsub IGHcoding_exons.sh 

bash /users/lcollado/R/x86_64-pc-linux-gnu-library/3.4.x/recount.bwtool/extdata/jhpce/run_rse.sh -r "/dcl01/lieber/ajaffe/lab/insp/IGHcoding_exons/IGHcoding_exons.Rdata" -s "sumsIGHcoding_exons" -c 1
