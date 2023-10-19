##############################################
# Metagenome Assembly Pipeline
##############################################

##### Merging paired-end reads with BBTools (bbmerge.sh) #####

module load iccifort/2020.1.217
module load BBMap/38.87

cd $SCRATCH/BBMap/bbmerge

while IFS= read -r i || [[ -n "$i" ]]
do 
bbmerge-auto.sh in1=$HOME/amrplusplus/NonHostReads/$i.non.host.R1.fastq.gz in2=$HOME/amrplusplus/NonHostReads/$i.non.host.R2.fastq.gz out=./$i.merged.fq outu=./$i.unmerged.fq ihist=./$i.ihist.txt ecct extend2=20 iterations=5
done < $HOME/ERIN_samples_IDs.txt;


##### Perform metagenome assembly with MEGAHit #####

module load MEGAHIT/1.2.9

cd $SCRATCH

while IFS= read -r i || [[ -n "$i" ]]
do 
megahit -1 $SCRATCH/ERIN_amrplusplus/NonHostReads/$i.non.host.R1.fastq.gz -2 $SCRATCH/ERIN_amrplusplus/NonHostReads/$i.non.host.R2.fastq.gz -r ./BBMap/bbmerge/$i.merged.fq --out-dir ./MEGAHIT/$i --out-prefix $i;
megahit_toolkit contig2fastg 99 ./MEGAHIT/$i/intermediate_contigs/k99.contigs.fa > ./MEGAHIT/$i/$i.k99.fastg 
done < $HOME/ERIN_samples_IDs_SUB2.txt;


##### Anvi'o Workflow - Assembled Contigs #####

module load iccifort/2019.5.281  impi/2018.5.288
module load anvio/7.0-Python-3.7.4

cd $SCRATCH/Anvio

while IFS= read -r i || [[ -n "$i" ]]
do 
# Reformat FASTA files (contigs) for use with Anvi'o
anvi-script-reformat-fasta ../MEGAHIT/$i.contigs.fa -o ./reformat/$i.contigs_fixed.fa -l 500 --simplify-names -r ./reformat/$i.report.txt

#Configure contig database with Anvi'o
anvi-gen-contigs-database -f ./reformat/$i.contigs_fixed.fa -o ./contigs_database/$i.contigs.db -n ./contigs_database/$i;

# Search for HMMs against the contig databases for each sample
anvi-run-hmms -c ./contigs_database/$i.contigs.db --num-threads 4;

# Export sequences of gene calls from the contig databases (after HMM discovery)
anvi-get-sequences-for-gene-calls -c ./contigs_database/$i.contigs.db --get-aa-sequences -o ./function/gene_calls/$i.aa_sequences.fa;
done < $HOME/ERIN_samples_IDs.txt;

# Note: Gene calls were used in the subsequent ACC analysis to identify ARGs on contigs


##### Assembly QC with MetaQUAST #####

cd $SCRATCH

module load QUAST/5.0.0

while IFS= read -r i || [[ -n "$i" ]]
do 
python /opt/software/QUAST/5.0.0/metaquast.py ./Anvio/reformat/$i.contigs_fixed.fa -f -o ./QUAST/$i/ -1 ./ERIN_amrplusplus/NonHostReads/$i.non.host.R1.fastq.gz -2 ./ERIN_amrplusplus/NonHostReads/$i.non.host.R2.fastq.gz 
done < $HOME/ERIN_samples_IDs.txt;












