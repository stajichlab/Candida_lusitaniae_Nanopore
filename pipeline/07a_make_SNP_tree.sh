#!/usr/bin/bash

#SBATCH --mem=24gb --ntasks 16 --nodes 1
#SBATCH --time=24:00:00
#SBATCH -J makeTree --out logs/make_tree.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [[ -f config.txt ]]; then
	source config.txt
else
	echo "Need a config.txt"
	exit
fi

if [[ -z $REFNAME ]]; then
	REFNAME=REF
fi
module load bcftools
module load IQ-TREE
module load fasttree
mkdir $TREEDIR
root=$FINALVCF/$PREFIX.selected.SNP
FAS=$TREEDIR/$PREFIX.fasaln
if [ -f $root.vcf ]; then
	module load tabix
	bgzip $root.vcf
fi
vcf=$root.vcf.gz
tab=$root.bcftools.tab
if [ ! -f $tab ]; then
	bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $tab
fi
if [ ! -f $FAS ]; then
    printf '>'$REFNAME'\n' > $FAS  
    bcftools query -f '%REF' ${vcf} >> $FAS
    printf '\n' >> $FAS

    for samp in $(bcftools query -l ${vcf} | grep -v -P '^CL_\d+'); do
	printf '>'${samp}'\n'
	bcftools query -s ${samp} -f '[%TGT]' ${vcf}
	printf '\n'
    done >> $FAS
fi
if [ ! -f $TREEDIR/$PREFIX.fasttree.tre ]; then
 FastTreeMP -gtr -gamma -nt < $FAS > $TREEDIR/$PREFIX.fasttree.tre
fi

iqtree-omp -nt $CPU -s $FAS -m GTR+ASC -b 100
