#!/usr/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 16G -p short
#SBATCH --job-name=GATK.select_filter
#SBATCH --output=logs/GATK.select_filter.%A.log

#Filters and selects for high quality variants, splitting them up into INDELONLY and SNPONLY files. Adjust filter parameters as needed

module load gatk/3.8

CONFIG=config.txt

if [[ -f $CONFIG ]]; then
	source $CONFIG
else
	echo "Expected a $CONFIG"
	exit
fi
INFILE=$FINALVCF/$PREFIX.all.vcf
SNP=$FINALVCF/$PREFIX.SNP.vcf
INDEL=$FINALVCF/$PREFIX.INDEL.vcf
FILTERSNP=$FINALVCF/$PREFIX.filtered.SNP.vcf
FILTERINDEL=$FINALVCF/$PREFIX.filtered.INDEL.vcf
FINALSNP=$FINALVCF/$PREFIX.selected.SNP.vcf
FINALINDEL=$FINALVCF/$PREFIX.selected.INDEL.vcf
GENOME=$GENOMEFOLDER/$GENOMENAME.fasta

## if [[ $GATK == "3" ]]; then
if [ ! -f $SNP ]; then
 java -Xmx3g -jar $GATK \
 -T SelectVariants \
 -R $GENOME \
 --variant $INFILE \
 -o $SNP \
 -env \
 -ef \
 -restrictAllelesTo BIALLELIC \
 -selectType SNP
fi

if [ ! -f $INDEL ]; then
 java -Xmx3g -jar $GATK \
 -T SelectVariants \
 -R $GENOME \
 --variant $INFILE \
 -o $INDEL \
 -env \
 -ef \
 -selectType INDEL -selectType MIXED -selectType MNP
fi

if [ ! -f $FILTERSNP ]; then
 java -Xmx3g -jar $GATK \
 -T VariantFiltration -o $FILTERSNP \
 --variant $SNP -R $GENOME \
 --clusterWindowSize 10  -filter "QD<2.0" -filterName QualByDepth \
 -filter "MQ<40.0" -filterName MapQual \
 -filter "QUAL<100" -filterName QScore \
 -filter "FS>60.0" -filterName FisherStrandBias \
 -filter "ReadPosRankSum<-8.0" -filterName ReadPosRank \
 --missingValuesInExpressionsShouldEvaluateAsFailing 

#-filter "HaplotypeScore > 13.0" -filterName HaplotypeScore
#-filter "MQ0>=10 && ((MQ0 / (1.0 * DP)) > 0.1)" -filterName MapQualRatio \
fi

if [ ! -f $FILTERINDEL ]; then
 java -Xmx3g -jar $GATK \
 -T VariantFiltration -o $FILTERINDEL \
 --variant $INDEL -R $GENOME \
 --clusterWindowSize 10  -filter "QD<2.0" -filterName QualByDepth \
 -filter "MQRankSum < -12.5" -filterName MapQualityRankSum \
 -filter "SOR > 4.0" -filterName StrandOddsRatio \
 -filter "FS>200.0" -filterName FisherStrandBias \
 -filter "InbreedingCoeff<-0.8" -filterName InbreedCoef \
 -filter "ReadPosRankSum<-20.0" -filterName ReadPosRank 
fi

if [ ! -f $FINALSNP ]; then
 java -Xmx16g -jar $GATK \
   -R $GENOME \
   -T SelectVariants \
   --variant $FILTERSNP \
   -o $FINALSNP \
   -env \
   -ef \
   --excludeFiltered
fi

if [ ! -f $FINALINDEL ]; then
 java -Xmx16g -jar $GATK \
   -R $GENOME \
   -T SelectVariants \
   --variant $FILTERINDEL \
   -o $FINALINDEL \
   --excludeFiltered 
fi

