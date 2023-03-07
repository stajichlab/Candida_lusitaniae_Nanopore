#!/usr/bin/bash -l
#SBATCH -p short -c 16 --mem 16gb -N 1 -n 1 --out logs/genespace_run.%A.log
module load orthofinder
module load MCScanX
module load R

Rscript genespace.R

