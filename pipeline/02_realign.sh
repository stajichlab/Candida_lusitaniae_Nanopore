#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH  --mem=32G
#SBATCH  --time=36:00:00
#SBATCH --job-name realign
#SBATCH --output=logs/realign.%a.log

module load java/8
module load gatk/3.7
module load picard

RGCENTER=MyCenter
RGPLATFORM=Illumina

CONFIG=config.txt

if [ -f $CONFIG ]; then
    source $CONFIG
fi

MEM=32g
GENOMEIDX=$GENOMEFOLDER/$GENOMENAME
KNOWNSITES=
if [ ! -f $GENOMEFOLDER/$GENOMENAME.dict ]; then
    picard CreateSequenceDictionary R=$GENOMEIDX O=$GENOMEFOLDER/$GENOMENAME.dict SPECIES=Candida_lusitaniae TRUNCATE_NAMES_AT_WHITESPACE=true
fi

if [ ! $CPU ]; then
    CPU=1
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	echo "Need a number via slurm array or on the cmdline"
	exit
    fi
fi

IFS=,
sed -n ${N}p $SAMPLESINFO | while read SAMPLE READ1 READ2 CTR TYPE
do
    if [ $TYPE == "Pool" ]; then
	echo "Skipping Pooled samples ($SAMPLE): Type=$TYPE"
	exit
    fi
    if [ ! -e $ALNFOLDER/$SAMPLE.DD.bam ]; then
	echo "Missing $ALNFOLDER/$SAMPLE.DD.bam - re-run step 1 with $N"
	exit
    fi
    if [ ! -e $ALNFOLDER/$SAMPLE.DD.bai ]; then
 	picard BuildBamIndex I=$ALNFOLDER/$SAMPLE.DD.bam TMP_DIR=/scratch
    fi

    if [ ! -e $ALNFOLDER/$SAMPLE.intervals ]; then 
 	java -Xmx${MEM} -jar $GATK \
   	    -T RealignerTargetCreator \
   	    -R $GENOMEIDX.fasta \
   	    -I $ALNFOLDER/$SAMPLE.DD.bam \
   	    -o $ALNFOLDER/$SAMPLE.intervals
    fi
    
    if [ ! -e $ALNFOLDER/$SAMPLE.realign.bam ]; then
	java -Xmx$MEM -jar $GATK \
   	    -T IndelRealigner \
   	    -R $GENOMEIDX.fasta \
   	    -I $ALNFOLDER/$SAMPLE.DD.bam \
   	    -targetIntervals $ALNFOLDER/$SAMPLE.intervals \
   	    -o $ALNFOLDER/$SAMPLE.realign.bam
    fi
    
    if [ ! -e $KNOWNSITES]; then
	if [ ! -f $ALNFOLDER/$SAMPLE.recal.grp ]; then
 	    java -Xmx$MEM -jar $GATK \
		-T BaseRecalibrator \
		-R $GENOMEIDX.fasta \
		-I $ALNFOLDER/$SAMPLE.realign.bam \
		--knownSites $KNOWNSITES \
		-o $ALNFOLDER/$SAMPLE.recal.grp
	fi
	if [ ! -f $ALNFOLDER/$SAMPLE.recal.bam ]; then
 	    java -Xmx$MEM -jar $GATK \
		-T PrintReads \
		-R $GENOMEIDX.fasta \
		-I $ALNFOLDER/$SAMPLE.realign.bam \
		-BQSR $ALNFOLDER/$SAMPLE.recal.grp \
		-o $ALNFOLDER/$SAMPLE.recal.bam
	fi
    fi
done
