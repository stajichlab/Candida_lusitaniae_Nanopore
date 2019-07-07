#!/usr/bin/bash
#SBATCH --mem 32G --ntasks 8 --nodes 1 -J ClusAln
#SBATCH --out logs/Clus.bwa.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard
module load gatk/3.7
module load samtools/1.9
module load bcftools/1.9
MEM=32g

TOPOUTDIR=tmp
ALNFOLDER=cram
HTCEXT=cram
HTCFORMAT=cram

if [ -f config.txt ]; then
    source config.txt
fi
if [ -z $REFGENOME ]; then
    echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
    exit
fi

if [ ! -f $REFGENOME.dict ]; then
    echo "NEED a $REFGENOME.dict - make sure 00_index.sh is run"
fi

mkdir -p $TOPOUTDIR
SAMPFILE=samples.csv
TEMP=/scratch

N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi


if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then 
 echo "need to provide a number by --array or cmdline"
 exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
echo "$N $MAX for $SAMPFILE"
if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
cat $SAMPFILE | sed -n ${N}p | while read SAMPLE PAIR1 PAIR2 CTR TYPE POPULATION
do
	echo "processing $SAMPLE $PAIR1 from file ($N)"
   # skip pool samples for now all together
   if [[ $TYPE != "Monoisolate" ]]; then
	   echo "skipping because type is $TYPE"
	   continue
   fi
   PAIR1=$FASTQFOLDER/$PAIR1
   PAIR2=$FASTQFOLDER/$PAIR2
  STRAIN=$(echo "$SAMPLE" | perl -p -e 's/ +/_/g; s/\//-/g;')
  SAMFILE=$TOPOUTDIR/$STRAIN.unsrt.sam
  SRTED=$TOPOUTDIR/${STRAIN}.srt.bam
  DDFILE=$TOPOUTDIR/${STRAIN}.DD.bam
  REALIGN=$TOPOUTDIR/${STRAIN}.realign.bam
  INTERVALS=$TOPOUTDIR/${STRAIN}.intervals
  FINALFILE=$ALNFOLDER/${STRAIN}.$HTCEXT    
  CENTER=$(echo $CTR | perl -p -e 's/ +/_/g')
  READGROUP="@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$STRAIN\tPL:illumina\tCN:$CENTER"
  echo "RG=$READGROUP"
  if [ ! -f $FINALFILE ]; then
      if [ ! -f $DDFILE ]; then
	  if [ ! -f $SRTED ]; then
	      if [ -e $PAIR2 ]; then
		  echo "RUNNING paired-end bwa"
		  bwa mem -t $CPU -R $READGROUP $REFGENOME $PAIR1 $PAIR2 > $SAMFILE
		  samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${STRAIN}.fixmate.bam
		  samtools sort --threads $CPU -O bam -o  $SRTED -T $TEMP $TEMP/${STRAIN}.fixmate.bam
		  /usr/bin/rm $TEMP/${STRAIN}.fixmate.bam $SAMFILE
	      elif [ -e $PAIR1 ]; then
		  echo "RUNNING unpaired bwa"
    		  bwa mem -t $CPU -R $READGROUP $REFGENOME $PAIR1 | samtools sort -@ $CPU -O bam -T $TEMP -o $SRTED
	      else
		  echo "NO $PAIR1 and no $PAIR2?"
		  exit
	      fi
	  fi # SRTED file exists or was created by this block
	  
	  echo "picard MarkDuplicates I=$SRTED O=$DDFILE \
	      METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null"

	  picard MarkDuplicates I=$SRTED O=$DDFILE \
	      METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null
	  if [ -f $DDFILE ]; then
	      rm -f $SRTED
	  fi
	  if [ ! -f $DDFILE.bai ]; then
	      picard BuildBamIndex I=$DDFILE TMP_DIR=/scratch
	  fi	    
      fi # DDFILE is created after this or already exists
      
      if [ ! -f $INTERVALS ]; then 
	  time java -Xmx$MEM -jar $GATK \
	      -T RealignerTargetCreator \
	      -R $REFGENOME \
	      -I $DDFILE \
	      -o $INTERVALS
      fi
      if [ ! -f $REALIGN ]; then
	  time java -Xmx$MEM -jar $GATK \
	      -T IndelRealigner \
	      -R $REFGENOME \
	      -I $DDFILE \
	      -targetIntervals $INTERVALS \
	      -o $REALIGN
      fi
      
      samtools view -O $HTCFORMAT --threads $CPU \
	  --reference $REFGENOME -o $FINALFILE $REALIGN
      samtools index $FINALFILE
      if [ -f $FINALFILE ]; then
	  rm -f $DDFILE $REALIGN
	  rm -f $(echo $REALIGN | sed 's/bam$/bai/')
	  rm -f $(echo $DDFILE | sed 's/bam$/bai/')
	  rm -f $INTERVALS
      fi
    else
	echo "$FINALFILE already exists"
    fi #FINALFILE created or already exists  
done
