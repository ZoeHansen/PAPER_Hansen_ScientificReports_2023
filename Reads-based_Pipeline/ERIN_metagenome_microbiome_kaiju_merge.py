# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 17:25:12 2021

@author: hanse
Last modified: 2021.10.25
"""

# This script is designed to compile all of the taxonomic classification
# completed by Kaiju into a single spreadsheet for analysis. 

# With the most recent run of Kaiju, this will need to be completed for 
# each separate taxonomic rank

#####################################################################

import pandas as pd
import os
import numpy as np

rootdir = r'D://HPCC/ERIN_Kaiju_taxonomic_classification_READS_102021/'
os.chdir(rootdir)
    
a=open('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/sampleIDs_all.txt')
a1 = a.read().splitlines()

# First we will address the 'kaiju_reads_summary_all.tsv' files.

for i in a1:
    kaiju_file = pd.read_csv(r"".join(i)+'.kaiju_SPECIES.tsv', sep='\t', header = None)
    kaiju_files = pd.DataFrame(kaiju_file)
    kaiju_files.columns = kaiju_files.iloc[0]
    kaiju_files.drop(kaiju_files.index[0], inplace=True)
    kaiju_files.columns = ['File','Percent',''+"".join(i)+'.reads', 'Taxon_ID','Taxon_name']
    kaiju_files.drop("File", axis=1, inplace=True) 
    kaiju_files.drop('Taxon_ID', axis=1, inplace=True)
    kaiju_files.to_csv('CSVs/'+''.join(i+'_kaiju_SPECIES.csv'), sep=',', index=False)

    print(i)  
  
a.close()

# Merge the classification data across all samples
# We can use the 'a1' variable to loop through the sample names

# Create a "starter" for the file merge
rootdir = r'D://HPCC/ERIN_Kaiju_taxonomic_classification_READS_102021/CSVs/'
os.chdir(rootdir)

merged = pd.DataFrame()
c = ['Percent','Taxon_name', 'ER0043','ER0073']
a2 = a1[2:]

main = pd.read_csv('ER0043_kaiju_SPECIES.csv', sep = ',', header = 0) 
main1=pd.DataFrame(main)
main1.drop('Percent', axis=1, inplace=True)
main1_summed=main1.groupby(by=['Taxon_name']).sum()
main1_summed.reset_index(inplace=True)
main1_summed.columns=['Taxon_name','ER0043']

main2 = pd.read_csv('ER0073_kaiju_SPECIES.csv', sep = ',', header = 0) 
main2=pd.DataFrame(main2)
main2.drop('Percent', axis=1, inplace=True)
main2_summed=main2.groupby(by=['Taxon_name']).sum()
main2_summed.reset_index(inplace=True)
main2_summed.columns=['Taxon_name','ER0073']

merged = pd.merge(main1_summed, main2_summed, on = 'Taxon_name', how = 'outer')

#Merge the rest of the samples:

for i in a2:
    main3 = pd.read_csv('' + ''.join(i) + '_kaiju_SPECIES.csv', sep = ',', header = 0)
    main3 = pd.DataFrame(main3)
    main3.drop('Percent', axis=1, inplace=True)
    main3_summed=main3.groupby(by=['Taxon_name']).sum()
    main3_summed.reset_index(inplace=True)
    main3_summed.columns=['Taxon_name',''+''.join(i)]
    
    merged = pd.merge(merged, main3_summed, on = 'Taxon_name', how = 'outer')
    merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(merged)

merged.to_csv('D://HPCC/ERIN_Kaiju_taxonomic_classification_READS_102021/ERIN_CaseControlFollow_kaiju_SPECIES.csv', sep = ',', index = False)


