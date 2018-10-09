#!/usr/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem=96G
#SBATCH --job-name=GATK.GVCFGeno
#SBATCH --output=logs/GATK.GVCFGeno.%A.log
#SBATCH --time=12:00:00

#Takes each individual sample vcf from Haplotype Caller step and combines it into single, combined vcf
MEM=96g #Requires large amount of memory. Adjust according to existing resources
module load picard
module load tabix
module load gatk/3.8
CONFIG=config.txt

if [ -f $CONFIG ]; then
    source $CONFIG
else
	echo "Expected a config file $CONFIG"
	exit
fi
GENOMEIDX=$GENOMEFOLDER/$GENOMENAME.fasta
KNOWNSITES=
OUT=$FINALVCF/$PREFIX.all.vcf
mkdir -p $FINALVCF
CPU=1

if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
if [[ $(ls $VARIANTFOLDER/*.g.vcf | wc -l | awk '{print $1}') -gt "0" ]]; then
	parallel -j 16 bgzip {} ::: $VARIANTFOLDER/*.g.vcf
	parallel -j 16 tabix -f {} ::: $VARIANTFOLDER/*.g.vcf.gz
fi

N=$(ls $VARIANTFOLDER/*.g.vcf.gz | sort | perl -p -e 's/\n/ /; s/(\S+)/-V $1/') #Lists each sample vcf by -V sample1.vcf -V sample2.vcf...

java -Xmx$MEM -jar $GATK \
    -T GenotypeGVCFs \
    -R $GENOMEIDX \
    $N \
    --max_alternate_alleles 3 \
    -o $OUT \
    -nt $CPU  
