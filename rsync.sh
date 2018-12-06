#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=800G
#$ -N rsync_insp
#$ -o ./rsync.txt
#$ -e ./rsync.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Copy with rsync
rsync -av /dcl01/lieber/ajaffe/lab/insp/ /dcl01/ajaffe/data/lab/insp/
## Copy faster (first time):
# https://superuser.com/questions/109780/how-to-speed-up-rsync
#cp -a /dcl01/lieber/ajaffe/lab/insp/ /dcl01/ajaffe/data/lab/

echo "**** Job ends ****"
date
