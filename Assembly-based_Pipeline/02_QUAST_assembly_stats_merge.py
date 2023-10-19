# -*- coding: utf-8 -*-
"""
@author: Zoe Hansen
Last Modified: 2022.01.11
"""

### This script will read in the output TSV files from QUAST, which holds information about the quality of 
### assembly by MEGAHIT. The output will then be merged for easier access/visualization

###########################
# Create CSV files for each sample 
###########################

#Read in TXT files to CSV files (try to remove columns, rename columns in the process)
import pandas as pd
import os
import numpy as np

rootdir = r'D://HPCC/QUAST_assembly_stats/'
os.chdir(rootdir)
    
a=open('samples_IDs_clean.txt')
a1 = a.read().splitlines()

for i in a1:
    file = pd.read_csv(r"".join(i)+'.report.tsv', sep='\t', header = 0)
    quast = pd.DataFrame(file)
    quast.columns = ['Assembly_stats', ''.join(i)]
    quast.to_csv('D://HPCC/QUAST_assembly_stats/CSV_files/'+"".join(i)+'.quast_report.csv', sep = ',', index=False)
    print(i)  
  
##################################
# Merge the CSV files  
##################################

# First, we will merge the 'ARGs_per_genera' files
dir1 = r'D://HPCC/QUAST_assembly_stats/CSV_files/'
os.chdir(dir1)

# Create a new variable 'a2' excluding our first two samples which initiate the merged file

a2 = a1[2:]

quast_main = pd.read_csv('ER0003.quast_report.csv', sep = ',', header = 0) 
quast_main1=pd.DataFrame(quast_main)

quast_main2 = pd.read_csv('ER0043.quast_report.csv', sep = ',', header = 0) 
quast_main2=pd.DataFrame(quast_main2)

quast_merged = pd.merge(quast_main1, quast_main2, on = ['Assembly_stats'], how = 'outer')

#Merge the rest of the samples:

for i in a2:
    quast_main3 = pd.read_csv(''.join(i) + '.quast_report.csv', sep = ',', header = 0)
    quast_main3 = pd.DataFrame(quast_main3)

    quast_merged = pd.merge(quast_merged, quast_main3, on = ['Assembly_stats'], how = 'outer')
    quast_merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(quast_merged)

quast_merged.to_csv('QUAST_assembly_statistics.csv', sep = ',', index = False)


