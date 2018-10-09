#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --time 2:00:00 -p short --mem 64G --out mosdepth.parallel.log
#SBATCH -J modepth
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
module unload python/2.7.5
mkdir -p coverage/mosdepth
export PATH="/bigdata/stajichlab/jstajich/miniconda3/bin:$PATH"

WINDOW=5000
parallel --jobs $CPU mosdepth -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:aln\/:coverage/mosdepth/:; s:\.realign\.bam:.${WINDOW}bp: =}" {} ::: aln/*.realign.bam

WINDOW=10000
parallel --jobs $CPU mosdepth -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:aln\/:coverage/mosdepth/:; s:\.realign\.bam:.${WINDOW}bp: =}" {} ::: aln/*.realign.bam

WINDOW=20000
parallel --jobs $CPU mosdepth -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:aln\/:coverage/mosdepth/:; s:\.realign\.bam:.${WINDOW}bp: =}" {} ::: aln/*.realign.bam

bash scripts/mosdepth_prep_ggplot.sh
