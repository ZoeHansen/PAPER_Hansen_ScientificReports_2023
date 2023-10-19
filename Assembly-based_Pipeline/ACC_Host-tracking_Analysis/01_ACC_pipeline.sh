############################################
# ARG-carrying Contig Analysis Pipeline
############################################

##### Annotation of ARGs on contigs using HMD-ARG and Diamond #####

module load DIAMOND/2.0.1

# Make the DIAMOND database
cd $HOME/Database/hmd-arg/
diamond makedb --in arg_v5.fasta -d hmd-arg

# Run the alignment
cd $HOME/ACC_pipeline/hmd-arg_alignment

while IFS= read -r i || [[ -n "$i" ]]
do 
diamond blastp \
      -d $HOME/Database/hmd-arg/hmd-arg.dmnd \
      -q /$SCRATCH/Anvio/function/gene_calls/$i.aa_sequences.fa.gz \
      -o ./$i.hmdarg_matches.sam \
      -f 101 \
      --evalue 0.00001 \
      --id 80 \
      -k 1 \
      --max-hsps 1 \
      --verbose \
      --threads 50
echo "HMD-ARG alignment done: $i"
done < $HOME/ERIN_samples_IDs.txt;


##### Identify which contigs contain ARGs #####

cd $HOME/ACC_pipeline

while IFS= read -r i || [[ -n "$i" ]]
do 
awk '/antibiotic/{print $1}' ./hmd-arg_alignment/$i.hmdarg_matches.sam > ./arg_contigs_all/$i.arg_contigs_all.lst
uniq -u ./arg_contigs_all/$i.arg_contigs_all.lst ./arg_contigs_unique/$i.arg_contigs_uniq.lst
echo "ACC extraction complete: $i"
done < $HOME/ERIN_samples_IDs.txt;


##### Use "seqtk" to extract all ARG-carrying contigs into a new FASTA file #####

module load GCC/9.3.0
module load seqtk/1.3

cd $SCRATCH/ACC_pipeline/

while IFS= read -r i || [[ -n "$i" ]]
do 
seqtk subseq $SCRATCH/Anvio/function/gene_calls/$i.aa_sequences.fa.gz ./arg_contigs_unique/$i.arg_contigs_uniq.lst > ./arg-carrying_contigs/$i.acc.fasta
echo "ACC extraction done: $i"
done < $HOME/ERIN_samples_IDs.txt;


##### Perform BLASTP annotation to identify taxa associated with ACCs #####

module load GCC/9.3.0  OpenMPI/4.0.3
module load BLAST+/2.10.1
export BLASTDB=/mnt/research/common-data/Bio/blastdb:/mnt/research/common-data/Bio/blastdb/v5:$BLASTDB

cd $SCRATCH/ACC_pipeline/

while IFS= read -r i || [[ -n "$i" ]]
do 
blastp \
   -query ./arg-carrying_contigs/$i.acc.fasta \
   -db nr \
   -out ./blastp_acc_output/$i.acc_blastp.txt \
   -evalue 0.00001 \
   -num_threads 50 \
   -max_target_seqs 50 \
   -outfmt "6 qseqid qgi qacc qlen sallseqid sallgi sacc slen qstart qend sstart send evalue bitscore score length mismatch gapopen pident nident btop staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp"
echo "blastp annotation done: $i"
done < $HOME/ERIN_samples_IDs_SUB.txt;

