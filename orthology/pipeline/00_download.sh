mkdir -p gff dna cds
curl -C- -o gff/Clus_L17.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/675/555/GCA_003675555.2_ASM367555v2/GCA_003675555.2_ASM367555v2_genomic.gtf.gz
curl -C- -o dna/Clus_L17.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/675/555/GCA_003675555.2_ASM367555v2/GCA_003675555.2_ASM367555v2_genomic.fna.gz
curl -C- -o cds/Clus_L17.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/675/555/GCA_003675555.2_ASM367555v2/GCA_003675555.2_ASM367555v2_cds_from_genomic.fna.gz

curl -C- -o cds/Clus_ATCC42720.fasta https://fungidb.org/common/downloads/release-62/ClusitaniaeATCC42720/fasta/data/FungiDB-62_ClusitaniaeATCC42720_AnnotatedCDSs.fasta
curl -C- -o dna/Clus_ATCC42720.fasta https://fungidb.org/common/downloads/release-62/ClusitaniaeATCC42720/fasta/data/FungiDB-62_ClusitaniaeATCC42720_Genome.fasta
curl -C- -o gff/Clus_ATCC42720.gff https://fungidb.org/common/downloads/release-62/ClusitaniaeATCC42720/gff/data/FungiDB-62_ClusitaniaeATCC42720.gff

# no alt splicing so we can take easy way out on this
zcat cds/Clus_L17.fasta.gz | perl -p -e 's/>(\S+)(.+)\[locus_tag=([^]]+)\]/>$3 $1$2/' > cds/Clus_L17.cds
perl -p -e 's/>(\S+)\-t\d+_\d+/>$1/' cds/Clus_ATCC42720.fasta > cds/Clus_ATCC42720.cds

mkdir peptide bed
module load emboss
transeq -table 12 -trim cds/Clus_ATCC42720.cds peptide/Clus_ATCC42720.fa
transeq -table 12 -trim cds/Clus_L17.cds peptide/Clus_L17.fa
perl -i -p -e 's/>(\S+)_1.+/>$1/' peptide/*.fa

grep -P "\tgene\t" gff/Clus_L17.gtf | cut -f 1,4,5,9 | grep -v 'gene_biotype "tRNA"' | perl -p -e 's/gene_id\s\"([^\"]+)\";.+/$1/' > bed/Clus_L17.bed
grep -P "\tprotein_coding_gene\t" gff/Clus_ATCC42720.gff | cut -f 1,4,5,9 | perl -p -e 's/ID=([^;]+);.+/$1/' > bed/Clus_ATCC42720.bed


