#!/usr/bin/bash

#SBATCH --mem=32gb --ntasks 2 --nodes 1
#SBATCH --time=24:00:00
#SBATCH -J fasalnTree


if [[ -f config.txt ]]; then
	source config.txt
else
	echo "Need a config.txt"
	exit
fi
module load vcftools
module load IQ-TREE
module load fastree

if [[ -z $REFNAME ]]; then
	REFNAME="REF"
fi

root=$FINALVCF/$PREFIX.selected.INDEL
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
    printf '>'$REF'\n' > $FAS
    bcftools query -f '%REF' ${vcf} >> $FAS
    printf '\n' >> $FAS

    for samp in $(bcftools query -l ${vcf} | grep -v -P '^CL_\d+'); do
	printf '>'${samp}'\n'
	bcftools query -s ${samp} -f '[%TGT]' ${vcf}
	printf '\n'
    done >> $FAS
fi

# INDEL tree
iqtree-omp -nt 2 -s indel_matrix.fasaln -m JC2 -b 100 -redo
