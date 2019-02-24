#!/bin/bash
#SBATCH --nodes 1 --ntasks 2 --mem 8G -p short --out logs/init.log

# EXPECTED VARIABLES
GENOMEFOLDER=genome
#GENOMEFASTA=C_lusitaniae_4606.v2.fasta
#GENOMENAME=ClusDS4606
CONFIG=config.txt

if [ -f $CONFIG ]; then
     source $CONFIG
fi

module load bwa/0.7.17
module load samtools/1.9
module load picard
mkdir -p logs

pushd $GENOMEFOLDER
FASTAFILE=$GENOMEFASTA
if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.bwt || $FASTAFILE -nt $FASTAFILE.bwt ]]; then
	bwa index $FASTAFILE
fi

DICT=$(basename $FASTAFILE .fasta)".dict"

if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	rm -f $DICT
	picard CreateSequenceDictionary R=$FASTAFILE O=$DICT
	ln -s $DICT $FASTAFILE.dict
fi

popd
