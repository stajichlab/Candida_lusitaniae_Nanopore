#!/usr/bin/bash 
#SBATCH --nodes 1 --ntasks 2 --mem 16G -p short --out logs/snpEff.log

module unload perl
module load perl/5.20.2
module load snpEff
module load bcftools/1.9
module load tabix

SNPEFFOUT=snpEff
SNPEFFGENOME=Clus_L1B
snpEffConfig=snpEff.config
GFFGENOME=Clavispora_lusitaniae_L1B.gff3
MEM=16g

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
	source config.txt
fi

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
	source config.txt
fi
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA
if [ -z $SNPEFFJAR ]; then
 echo "need to defined \$SNPEFFJAR in module or config.txt"
 exit
fi
if [ -z $SNPEFFDIR ]; then
 echo "need to defined \$SNPEFFDIR in module or config.txt"
 exit
fi
# could make this a confi

if [ -z $FINALVCF ]; then
	echo "need a FINALVCF variable in config.txt"
	exit
fi

mkdir -p $SNPEFFOUT
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
	rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
	echo "# $SNPEFFGENOME" >> $SNPEFFOUT/$snpEffConfig
  	echo "$SNPEFFGENOME.genome : Clavispora lusitaniae L1B" >> $SNPEFFOUT/$snpEffConfig
	chroms=$(awk '{print $1}' genome/Clavispora_lusitaniae_L1B.gff3 | grep -v '#' | sort | uniq | perl -p -e 's/\n/, /' | perl -p -e 's/,\s+$/\n/')
	
	echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig
	for chr in $(awk '{print $1}' genome/Clavispora_lusitaniae_L1B.gff3 | grep -v '#' | sort | uniq)
	do
	    echo -e "\t$SNPEFFGENOME.$chr.codonTable : codon.Alternative_Yeast_Nuclear"
	done

	#echo -e "\t$SNPEFFGENOME.mito_A_fumigatus_Af293.codonTable : Mold_Mitochondrial" >> $SNPEFFOUT/$snpEffConfig
	mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
	gzip -c $GFFGENOMEFILE > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
	rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa

	java -Xmx$MEM -jar $SNPEFFJAR build -datadir `pwd`/$SNPEFFOUT/data -c $SNPEFFOUT/$snpEffConfig -gff3 -v $SNPEFFGENOME
fi
pushd $SNPEFFOUT
COMBVCF="../$FINALVCF/$PREFIX.selected.SNP.vcf.gz ../$FINALVCF/$PREFIX.selected.INDEL.vcf.gz"

for n in $COMBVCF
do
 st=$(echo $n | perl -p -e 's/\.gz//')
 if [ ! -f $n ]; then
	 bgzip $st
	 tabix $n
 fi
done
INVCF=$PREFIX.combined_selected.SNP.vcf
OUTVCF=$PREFIX.snpEff.vcf
OUTTAB=$PREFIX.snpEff.tab
bcftools concat -a -d both -o $INVCF -O v $COMBVCF
java -Xmx$MEM -jar $SNPEFFJAR eff -dataDir `pwd`/data -v $SNPEFFGENOME $INVCF > $OUTVCF

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF > $OUTTAB
