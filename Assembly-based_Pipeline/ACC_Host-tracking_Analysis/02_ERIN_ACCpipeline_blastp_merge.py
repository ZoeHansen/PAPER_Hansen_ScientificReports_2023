# -*- coding: utf-8 -*-
"""
@author: Zoe Hansen
Last Modified: 2021.08.16
"""
### This script will first create CSV files and modify the blastp output from the ARG-carrying contigs pipeline.
### We will then identify the most-likely taxa identified in relation to ARGs on each contig within each sample.
### For this, a new CSV file will be produced for each sample with one taxa-ARG connection per contig.
### Upon completion of this ranking, we can merge the final files for downstream analysis.

###########################
# Create CSV files for each sample directly from BLASTP output
###########################

#Read in TXT files to CSV files (try to remove columns, rename columns in the process)
import pandas as pd
import os
import numpy as np

rootdir = r'D://HPCC/ACC_pipeline/BLASTP_annotations/TXT_files/'
os.chdir(rootdir)
    
a=open('ERIN_samples_IDs_clean.txt')
a1 = a.read().splitlines()

for i in a1:
    file = pd.read_csv(r"".join(i)+'.acc_blastp.txt', sep='\t', header = None)
    accs = pd.DataFrame(file)
    accs.columns = ['qseqid', 'qgi', 'qacc', 'qlen', 'sallseqid', 'sallgi', 'sacc', 'slen', 'qstart', 'qend', 'sstart', 
                    'send', 'evalue', 'bitscore', 'score', 'length', 'mismatch', 'gapopen', 'pident', 'nident', 
                    'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'qcovs',
                    'qcovhsp']
    accs.to_csv('D://HPCC/ACC_pipeline/BLASTP_annotations/CSV_files/'+"".join(i)+'.acc_blastp.csv', sep = ',', index=False)
    print(i)  
  
a.close()

##################################################################################

# The following code will be performed for a SINGLE sample to explain which data are
# being isolated in each step. The remaining samples will be considered/merged in a 
# for-loop below.

##################################
# Quantify the number of ACCs per sample
##################################

# This code will collapse all unique contig names for each sample and retrieve a count of 
# the total number of ACCs per sample. These counts will be concatenated into a merged file

newdir = r'D://HPCC/ACC_pipeline/BLASTP_annotations/CSV_files/'
os.chdir(newdir)
  
acc_csv = pd.read_csv(r'ER0003.acc_blastp.csv', sep=',', header = 0)
acc_csv1 = pd.DataFrame(acc_csv)
acc_quant = acc_csv1[['qseqid']]
acc_quant1 = acc_quant.groupby('qseqid').agg(''.join)
acc_quant1.reset_index(level=0, inplace=True)
acc_quant1.columns=['ER0003']
acc_quant2= pd.DataFrame(acc_quant1.count())
acc_quant2.reset_index(level=0, inplace=True)
acc_quant2.columns=['ER_ID','num_accs']

acc_csv2 = pd.read_csv(r'ER0043.acc_blastp.csv', sep=',', header = 0)
acc_csv3 = pd.DataFrame(acc_csv2)
acc_quant3 = acc_csv3[['qseqid']]
acc_quant4 = acc_quant3.groupby('qseqid').agg(''.join)
acc_quant4.reset_index(level=0, inplace=True)
acc_quant4.columns=['ER0043']
acc_quant5= pd.DataFrame(acc_quant4.count())
acc_quant5.reset_index(level=0, inplace=True)
acc_quant5.columns=['ER_ID','num_accs']

acc_quant_merged=pd.merge(acc_quant2, acc_quant5, how='outer', on=['ER_ID','num_accs'])

a2 = a1[2:]

for i in a2:
    acc_csv4 = pd.read_csv(r''.join(i)+'.acc_blastp.csv', sep=',', header = 0)
    acc_csv5 = pd.DataFrame(acc_csv4)
    acc_quant6 = acc_csv5[['qseqid']]
    acc_quant7 = acc_quant6.groupby('qseqid').agg(''.join)
    acc_quant7.reset_index(level=0, inplace=True)
    acc_quant7.columns=[''.join(i)]
    acc_quant8= pd.DataFrame(acc_quant7.count())
    acc_quant8.reset_index(level=0, inplace=True)
    acc_quant8.columns=['ER_ID','num_accs']
    acc_quant_merged=pd.merge(acc_quant_merged, acc_quant8, how='outer', on=['ER_ID','num_accs'])

print(acc_quant_merged)

acc_quant_merged.to_csv(r'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ACC_analysis/ERIN_total_ACCs_per_sample.csv', sep = ',', index = False)


##################################
# Rank the identification of various taxa on contigs
##################################

# The following code was used for a single sample; a larger code block below
# contains a for-loop which completes the following steps for all samples. 

# Note: I wanted to retain the single-sample code in case troubelshooting is needed

newdir = r'D://HPCC/ACC_pipeline/BLASTP_annotations/CSV_files/'
os.chdir(newdir)

blast = pd.read_csv('ER0003.acc_blastp.csv', sep=',', header=0)
blast1 = pd.DataFrame(blast)

# Remove unneeded columns
blast1.drop(blast1.loc[:,'qgi':'staxids'].columns, inplace=True, axis=1)
blast1.drop(blast1.loc[:,'scomnames':'sblastnames'].columns, inplace=True, axis=1)
blast1.drop(blast1.loc[:,'qcovs':'qcovhsp'].columns, inplace=True, axis=1)

# Cleanup our regex in the remaining columns
blast1.replace(';',' ', regex=True, inplace=True)

# Isolate the genus name in our blast hits
#blast1['genus']=blast1['sscinames']
blast1['genus'] = blast1['sscinames'].str.split().str.get(0)

# Remove the "sscinames" column since we extracted the genus information
blast1.drop('sscinames', inplace=True, axis=1)

# Find the cumulative sum of each genus per contig 
blast2=blast1.groupby('qseqid')['genus'].value_counts()

# Need to reset the index to obtain 3 columns: genus, contig ID, and number of hits
blast3=blast2.to_frame()
blast3.reset_index(level=0, inplace=True)
blast3.columns = ['contig_id','hits']
blast3.reset_index(inplace=True)

# Find the percentage per contig occupied by each genus:
blast2_pcts=blast2.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))

# Reset indices similar to above
blast3_pcts=blast2_pcts.to_frame()
blast3_pcts.reset_index(level=0, inplace=True)
blast3_pcts.columns = ['contig_id','pct']
blast3_pcts.reset_index(inplace=True)

# Merge these two dataframes to append percentages per contig assigned to each genus
blast_genus=pd.merge(blast3, blast3_pcts, how='outer', on=['contig_id','genus'])
blast_genus = blast_genus[['contig_id','genus','hits','pct']]


# Cleanup the ARG annotation and find the cumulative sum of each ARG per contig
blast1['arg']=blast1['stitle'].str.split('[').str[0]
blast1.drop('stitle', inplace=True, axis=1)
blast4=blast1.groupby('qseqid')['arg'].value_counts()

blast5=blast4.to_frame()
blast5.reset_index(level=0, inplace=True)
blast5.columns=['contig_id','hits']
blast5.reset_index(inplace=True)

blast4_pcts=blast4.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))

blast5_pcts=blast4_pcts.to_frame()
blast5_pcts.reset_index(level=0,inplace=True)
blast5_pcts.columns=['contig_id','pct']
blast5_pcts.reset_index(inplace=True)

blast_args=pd.merge(blast5, blast5_pcts, how='outer', on=['contig_id','arg'])
blast_args=blast_args[['contig_id','arg','hits','pct']]

### We now have two dataframes which contain the number of hits and percentages assigned to each taxa
### and ARG identified on our samples' contigs. 

blast_merged=pd.merge(blast_genus, blast_args, how='outer', on='contig_id')

# This output contains all ACCs, the number of hits/percent assigned to every genus and the number of hits/percent
# assigned to every ARG. This information may be useful to come back to if needed as a reference/to de-bug 
# something, but otherwise will be difficult to glean information from. 

# However, CSV files can be saved as a "checkpoint" in this code and for easier viewing/explanation.

# Write these dataframes to CSVs for each sample
blast_merged.to_csv('D://HPCC/ER0003_test_blastmerged.csv', sep=',', index=False)

#####################################
# Extracting Genera and ARGs of highest proportion per contig
#####################################

# Now, we will extract the genus or ARG with the greatest percentage on each contig
# Note: this will be done separately for genera and ARGs, so we will use the blast_genus and 
# blast_args dataframes sequentially

top_genus = pd.DataFrame(blast_genus)
top_genus=top_genus.loc[top_genus.groupby(['contig_id'])['pct'].idxmax()]
top_genus.reset_index(drop=True, inplace=True)
top_genus.columns=['contig_id','genus','genus_hits','genus_pct']

top_args = pd.DataFrame(blast_args)
top_args = top_args.loc[top_args.groupby(['contig_id'])['pct'].idxmax()]
top_args.reset_index(drop=True, inplace=True)
top_args.columns=(['contig_id','arg','arg_hits','arg_pct'])

top_merged = pd.merge(top_genus, top_args, how='outer', on='contig_id')

# In this output, we maintain the individual contig IDs, their top taxonomic assignment, and their top ARG
# assignment. This output will be written to separate CSV files to conserve this information, but is otherwise
# difficult to use in analysis. 

# The next section of this code is designed to further sort/manipulate these data to extract useful information
# that can be concatenated across all samples

# Write to CSVs
top_merged.to_csv('D://HPCC/ER0003_test_topmerge.csv', sep = ',', index = False)

##################################
# Sort by genus to view ARGs attributed to each genus
##################################

# Group-by genus for easier viewing of genus-ARG connections
hosts = top_merged.sort_values(by=['genus'])
hosts.reset_index(drop=True, inplace=True)
hosts2 = hosts.groupby('genus')['arg'].value_counts()

hosts3=hosts2.to_frame()
hosts3.reset_index(level=0, inplace=True)
hosts3.columns=['genus','hits']
hosts3.reset_index(level=0, inplace=True)
hosts3 = hosts3[['genus','arg','hits']]

# Determine percentages of ARGs assigned within each genus:
hosts2_pcts=hosts2.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))

hosts3_pcts=hosts2_pcts.to_frame()
hosts3_pcts.columns=['arg_pct']
hosts3_pcts.reset_index(level=0, inplace=True)
hosts3_pcts.columns=['genus','arg_pct']
hosts3_pcts.reset_index(level=0, inplace=True)
hosts3_pcts = hosts3_pcts[['genus','arg','arg_pct']]

host_merged=pd.merge(hosts3, hosts3_pcts, how='outer', on=['genus','arg'])

# In this output, the 'arg_pct' column shows the percentage of all ARGs assigned within a genus that is
# comprised of that particular ARG. 
#   For example, the MFS transporter registers a percentage of 8.7% in the Citrobacter genus within sample ER0003. 
#   This indicates that of all ARGs "assigned" to Citrobacter in ER0003, MFS transporters make up 8.7% of these.

# Write to CSV:
host_merged.to_csv('D://HPCC/ER0003_test_hostmerged.csv', sep=',', index=False)

##################################
# Find which genera comprise ACCs for each sample (percentages of total)
##################################

# Remove ARG information and sort by 'genus' column; obtain a "count" of these columns and divide by total 
# number of rows to achieve percentage of each sample's ACCs attributed to each genus

genera = pd.DataFrame(top_merged['genus'].value_counts())
genera.reset_index(level=0, inplace=True)
genera.columns = ['genus','hits']

# Determine percentage of ACCs attributed to each genus

genera['genus_percent'] = (genera['hits'] / genera['hits'].sum()) * 100 

# In this output, the 'genus_percent' column indicates the percentage of all genera assigned to ACCs in a sample. 
#   For example, Citrobacter registered a percentage of 39.58% in sample ER0003. This indicates that Citrobacter 
#   accounted for ~40% of all taxonomically-assigned ACCs in ER0003.

# Write to CSV:
genera.to_csv('D://HPCC/ER0003_test_generapct.csv', sep=',', index=False)

##########################################################################
##########################################################################

# For-loop that performs code above for all samples 

##########################################################################
##########################################################################

# Redefine 'a1' just in case it was lost prior    
a=open('ERIN_samples_IDs_clean.txt')
a1 = a.read().splitlines()

# Change directory to location of BLASTP CSV output files
newdir = r'D://HPCC/ACC_pipeline/BLASTP_annotations/CSV_files/'
os.chdir(newdir)

for i in a1:
    blast = pd.read_csv(r''.join(i)+'.acc_blastp.csv', sep=',', header=0)
    blast1 = pd.DataFrame(blast)
    blast1.drop(blast1.loc[:,'qgi':'staxids'].columns, inplace=True, axis=1)
    blast1.drop(blast1.loc[:,'scomnames':'sblastnames'].columns, inplace=True, axis=1)
    blast1.drop(blast1.loc[:,'qcovs':'qcovhsp'].columns, inplace=True, axis=1)
    blast1.replace(';',' ', regex=True, inplace=True)
    # GENUS
    blast1['genus'] = blast1['sscinames'].str.split().str.get(0)
    blast1.drop('sscinames', inplace=True, axis=1)
    blast2=blast1.groupby('qseqid')['genus'].value_counts()
    blast3=blast2.to_frame()
    blast3.reset_index(level=0, inplace=True)
    blast3.columns = ['contig_id','genus_hits']
    blast3.reset_index(inplace=True)
    blast2_pcts=blast2.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))
    blast3_pcts=blast2_pcts.to_frame()
    blast3_pcts.reset_index(level=0, inplace=True)
    blast3_pcts.columns = ['contig_id','genus_pct']
    blast3_pcts.reset_index(inplace=True)
    blast_genus=pd.merge(blast3, blast3_pcts, how='outer', on=['contig_id','genus'])
    blast_genus = blast_genus[['contig_id','genus','genus_hits','genus_pct']]
    # ARGs
    blast1['arg']=blast1['stitle'].str.split('[').str[0]
    blast1.drop('stitle', inplace=True, axis=1)
    blast4=blast1.groupby('qseqid')['arg'].value_counts()
    blast5=blast4.to_frame()
    blast5.reset_index(level=0, inplace=True)
    blast5.columns=['contig_id','arg_hits']
    blast5.reset_index(inplace=True)
    blast4_pcts=blast4.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))
    blast5_pcts=blast4_pcts.to_frame()
    blast5_pcts.reset_index(level=0,inplace=True)
    blast5_pcts.columns=['contig_id','arg_pct']
    blast5_pcts.reset_index(inplace=True)
    blast_args=pd.merge(blast5, blast5_pcts, how='outer', on=['contig_id','arg'])
    blast_args=blast_args[['contig_id','arg','arg_hits','arg_pct']]
    # MERGED
    blast_merged=pd.merge(blast_genus, blast_args, how='outer', on='contig_id')
    blast_merged.to_csv('D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_genera_arg_total_merged/'+"".join(i)+'_genera_arg_total_merged.csv', sep=',', index=False)
    print(i,': blast_merged')
    # TOP GENERA & ARGs
    top_genus = pd.DataFrame(blast_genus)
    top_genus=top_genus.loc[top_genus.groupby(['contig_id'])['genus_pct'].idxmax()]
    top_genus.reset_index(drop=True, inplace=True)
    top_genus.columns=['contig_id','genus','genus_hits','genus_pct']
    top_args = pd.DataFrame(blast_args)
    top_args = top_args.loc[top_args.groupby(['contig_id'])['arg_pct'].idxmax()]
    top_args.reset_index(drop=True, inplace=True)
    top_args.columns=(['contig_id','arg','arg_hits','arg_pct'])
    top_merged = pd.merge(top_genus, top_args, how='outer', on='contig_id')
    top_merged.to_csv('D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_top_genera_arg_merged/'+"".join(i)+'_top_genera_args.csv', sep = ',', index = False)
    print(i,': top_merged')
    # ARGs PER GENUS
    hosts = top_merged.sort_values(by=['genus'])
    hosts.reset_index(drop=True, inplace=True)
    hosts2 = hosts.groupby('genus')['arg'].value_counts()
    hosts3=hosts2.to_frame()
    hosts3.reset_index(level=0, inplace=True)
    hosts3.columns=['genus','hits']
    hosts3.reset_index(level=0, inplace=True)
    hosts3 = hosts3[['genus','arg','hits']]
    hosts2_pcts=hosts2.groupby(level=0).apply(lambda x: 100*x/float(x.sum()))
    hosts3_pcts=hosts2_pcts.to_frame()
    hosts3_pcts.columns=['arg_pct']
    hosts3_pcts.reset_index(level=0, inplace=True)
    hosts3_pcts.columns=['genus','arg_pct']
    hosts3_pcts.reset_index(level=0, inplace=True)
    hosts3_pcts = hosts3_pcts[['genus','arg','arg_pct']]
    host_merged=pd.merge(hosts3, hosts3_pcts, how='outer', on=['genus','arg'])
    host_merged.to_csv('D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_args_per_genus/'+"".join(i)+'_args_per_genus.csv', sep=',', index=False)
    print(i,': ARGs_per_genus')
    # GENERA PER SAMPLE
    genera = pd.DataFrame(top_merged['genus'].value_counts())
    genera.reset_index(level=0, inplace=True)
    genera.columns = ['genus','hits']
    genera['genus_percent'] = (genera['hits'] / genera['hits'].sum()) * 100 
    genera.to_csv('D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_genera_per_sample/'+"".join(i)+'_genera_per_sample.csv', sep=',', index=False)
    print(i,': genera_per_sample')
    
a.close()    

# At the end of this for loop, we should have 4 separate CSV files for each sample (designated by each checkpoint
# in this for-loop) 

##################################
# Merge the CSV files for "ARGs_per_genera" and "Genera_per_sample" files to obatin a comprehensive list 
##################################

# First, we will merge the 'ARGs_per_genera' files
dir1 = r'D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_args_per_genus/'
os.chdir(dir1)

# Create a new variable 'a2' excluding our first two samples which initiate the merged file

a2 = a1[2:]

apg_main = pd.read_csv('ER0003_args_per_genus.csv', sep = ',', header = 0) 
apg_main1=pd.DataFrame(apg_main)
apg_main1.drop('hits', axis=1, inplace=True)
apg_main1.columns = ['genus','arg','ER0003_pct']
#apg_main1.drop('arg_pct', axis = 1, inplace = True)
#apg_main1.columns = ['genus','arg','ER0003_hits']


apg_main2 = pd.read_csv('ER0043_args_per_genus.csv', sep = ',', header = 0) 
apg_main2=pd.DataFrame(apg_main2)
apg_main2.drop('hits', axis=1, inplace=True)
apg_main2.columns = ['genus','arg','ER00043_pct']
#apg_main2.drop('arg_pct', axis=1, inplace=True)
#apg_main2.columns = ['genus','arg','ER0043_hits']

apg_merged = pd.merge(apg_main1, apg_main2, on = ['genus','arg'], how = 'outer')

apg_merged.to_csv('D://HPCC/apg_merge_test.csv', sep=',', index=False)

#Merge the rest of the samples:

for i in a2:
    apg_main3 = pd.read_csv(''.join(i) + '_args_per_genus.csv', sep = ',', header = 0)
    apg_main3 = pd.DataFrame(apg_main3)
    apg_main3.drop('hits', axis=1, inplace=True)
    apg_main3.columns = ['genus','arg',''.join(i)+'_pct']
#    apg_main3.drop('arg_pct', axis=1, inplace=True)
#    apg_main3.columns = ['genus','arg',''.join(i)+'_hits']

    apg_merged = pd.merge(apg_merged, apg_main3, on = ['genus','arg'], how = 'outer')
    apg_merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(apg_merged)

apg_merged=apg_merged.sort_values(by=['genus','arg'])

apg_merged.to_csv('ERIN_args_per_genus_merged_pcts.csv', sep = ',', index = False)


###################################

# Now we will do the 'genera_per_sample' files

# First, we will merge the 'ARGs_per_genera' files
dir2 = r'D://HPCC/ACC_Pipeline/BLASTP_sorted_csvs/blastp_genera_per_sample/'
os.chdir(dir2)

# Create a new variable 'a2' excluding our first two samples which initiate the merged file

a2 = a1[2:]

gps_main = pd.read_csv('ER0003_genera_per_sample.csv', sep = ',', header = 0) 
gps_main1=pd.DataFrame(gps_main)
#gps_main1.drop('hits', axis=1, inplace=True)
#gps_main1.columns = ['genus','ER0003_pct']
gps_main1.drop('genus_percent', axis = 1, inplace = True)
gps_main1.columns = ['genus','ER0003_hits']


gps_main2 = pd.read_csv('ER0043_genera_per_sample.csv', sep = ',', header = 0) 
gps_main2=pd.DataFrame(gps_main2)
#gps_main2.drop('hits', axis=1, inplace=True)
#gps_main2.columns = ['genus','ER00043_pct']
gps_main2.drop('genus_percent', axis=1, inplace=True)
gps_main2.columns = ['genus','ER0043_hits']

gps_merged = pd.merge(gps_main1, gps_main2, on = ['genus'], how = 'outer')

#Merge the rest of the samples:

for i in a2:
    gps_main3 = pd.read_csv(''.join(i) + '_genera_per_sample.csv', sep = ',', header = 0)
    gps_main3 = pd.DataFrame(gps_main3)
#    gps_main3.drop('hits', axis=1, inplace=True)
#    gps_main3.columns = ['genus',''.join(i)+'_pct']
    gps_main3.drop('genus_percent', axis=1, inplace=True)
    gps_main3.columns = ['genus',''.join(i)+'_hits']

    gps_merged = pd.merge(gps_merged, gps_main3, on = ['genus'], how = 'outer')
    gps_merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(gps_merged)

gps_merged=gps_merged.sort_values(by=['genus'])

gps_merged.to_csv('ERIN_genera_per_sample_merged_hits.csv', sep = ',', index = False)


####################################################
# Investigating Percent Identity
####################################################

newdir = r'D://HPCC/ACC_pipeline/BLASTP_annotations/CSV_files/'
os.chdir(newdir)

acc_csv = pd.read_csv(r'ER0003.acc_blastp.csv', sep=',', header = 0)
acc_csv1 = pd.DataFrame(acc_csv)

acc_pctid = acc_csv1.groupby('qseqid')['pident'].mean().reset_index()
acc_pctid.insert(0,'contigs', range(1, 1+len(acc_pctid)))
acc_pctid.drop('qseqid', inplace=True, axis=1)
acc_pctid.columns = ['contigs','ER0003']


acc_csv2 = pd.read_csv(r'ER0043.acc_blastp.csv', sep=',', header = 0)
acc_csv3 = pd.DataFrame(acc_csv2)

acc_pctid2 = acc_csv3.groupby('qseqid')['pident'].mean().reset_index()
acc_pctid2.insert(0,'contigs', range(1, 1+len(acc_pctid2)))
acc_pctid2.drop('qseqid', inplace=True, axis=1)
acc_pctid2.columns = ['contigs','ER0043']


acc_pctid_merged=pd.merge(acc_pctid, acc_pctid2, how='outer', on=['contigs'])

a2 = a1[2:]

for i in a2:
    acc_csv4 = pd.read_csv(r''.join(i)+'.acc_blastp.csv', sep=',', header = 0)
    acc_csv5 = pd.DataFrame(acc_csv4)
    acc_pctid3 = acc_csv5.groupby('qseqid')['pident'].mean().reset_index()
    acc_pctid3.insert(0,'contigs', range(1, 1+len(acc_pctid3)))
    acc_pctid3.drop('qseqid', inplace=True, axis=1)
    acc_pctid3.columns=['contigs',''.join(i)]
    acc_pctid_merged=pd.merge(acc_pctid_merged, acc_pctid3, how='outer', on=['contigs'])

print(acc_pctid_merged)

acc_pctid_merged.to_csv(r'D://ACC_analysis/ERIN_avg_percent_id_per_contig.csv', sep = ',', index = False)




