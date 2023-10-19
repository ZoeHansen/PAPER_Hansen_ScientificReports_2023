################################################
# Kaiju Taxonomic Classification
################################################

### Create Kaiju database ("nr_euk" is the NCBI "nr" database including all proteins belonging to Archaea, bacteria, viruses, and fungi and other microbial eukaryotes

cd $HOME/Database/Kaiju

$HOME/kaiju/bin/kaiju-makedb -s nr_euk -t 100
wget http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2020-05-25.tgz

### Complete taxonomic classification with Kaiju

KAIJU_NR_DB=/mnt/ufs18/rs-033/common-data/kaiju_db_nr_euk/kaiju_db_nr_euk.fmi
KAIJU_NAMES=/mnt/ufs18/rs-033/common-data/kaiju_db_nr_euk/names.dmp
KAIJU_NODES=/mnt/ufs18/rs-033/common-data/kaiju_db_nr_euk/nodes.dmp

### Paired-End Read Annotation #####

cd /mnt/gs18/scratch/users/hansenzo/Kaiju/

while IFS= read -r i || [[ -n "$i" ]]
do 
$HOME/kaiju/bin/kaiju -t $KAIJU_NODES -f $KAIJU_NR_DB -i $HOME/amrplusplus/NonHostReads/$i.non.host.R1.fastq.gz -j $HOME/amrplusplus/NonHostReads/$i.non.host.R2.fastq.gz -o ./output/$i.kaiju.out -z 10 
done < $HOME/ERIN_samples_IDs.txt

while IFS= read -r i || [[ -n "$i" ]]
do 
$HOME/kaiju/bin/kaiju-addTaxonNames -t $KAIJU_NODES -n $KAIJU_NAMES -i ./output/$i.kaiju.out -o ./output/$i.kaiju.out.names -r superkingdom,phylum,class,order,family,genus,species -v;
done < $HOME/ERIN_samples_IDs.txt


### Create output tables at taxonomic levels of interest (Species, Genus, and Phylum here):

while IFS= read -r i || [[ -n "$i" ]]
do 
$HOME/kaiju/bin/kaiju2table -t $KAIJU_NODES -n $KAIJU_NAMES -o ./output/$i.kaiju_SPECIES.tsv ./output/$i.kaiju.out -r species -l superkingdom,phylum,class,order,family,genus,species -e -v;
$HOME/kaiju/bin/kaiju2table -t $KAIJU_NODES -n $KAIJU_NAMES -o ./output/$i.kaiju_GENUS.tsv ./output/$i.kaiju.out -r genus -l superkingdom,phylum,class,order,family,genus -e -v;
$HOME/kaiju/bin/kaiju2table -t $KAIJU_NODES -n $KAIJU_NAMES -o ./output/$i.kaiju_PHYLUM.tsv ./output/$i.kaiju.out -r phylum -l superkingdom,phylum -e -v;
done < $HOME/ERIN_samples_IDs.txt
