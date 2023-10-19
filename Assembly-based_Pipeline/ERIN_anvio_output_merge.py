# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 12:45:07 2021

@author: Zoe Hansen
Last Modified: 2021.08.11
"""

### This script is designed to merge the functional output from the Anvi'o pipeline.
### Specifically, this script will merge TXT files output from NCBI COG and/or
### SCG taxonomy annotations. 

#################################################

# NCBI COGS OUTPUT 

#################################################

#Read in .txt files to .csv files (try to remove columns, rename columns in the process)
import pandas as pd
import os
import numpy as np

rootdir = r'D://HPCC/Spring2021_Aim2_Pipeline/Anvio/06_2021_function_NCBI_COGS_output/'
os.chdir(rootdir)
    
a=open('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ERIN_samples_IDs_clean.txt')
a1 = a.read().splitlines()

for i in a1:
    file = pd.read_csv(r"".join(i)+'.ncbi_function.txt', sep='\t', header = None)
    funs = pd.DataFrame(file)
    funs.columns = funs.iloc[0]
    funs.drop(funs.index[0], inplace=True)
    funs.columns = ['gene_callers_id',	'source','accession','function',''.join(i)]
    funs2 = funs[funs.source == 'COG20_CATEGORY']
    funs2["".join(i)] = pd.to_numeric(funs2["".join(i)])
    funs2.sort_values(by=''.join(i), ascending=True, inplace=True)
    funs2.drop_duplicates(subset ="function", keep = 'first', inplace = True)
    funs2.to_csv('D://HPCC/Spring2021_Aim2_Pipeline/Anvio/06_2021_function_NCBI_COGS_output/'+"".join(i)+'_COG_categories.csv', sep = ',', index=False)
    print(i)  
  
a.close()


# Merge the resistome gene data (can use this code as a template for the mechanism, group, class data)
# We can use the 'a1' variable to loop through the sample names

# Create a "starter" for the file merge
merged = pd.DataFrame()
c = ['Function', 'ER0003','ER0043']
a2 = a1[2:]

main = pd.read_csv('ER0003_COG_categories.csv', sep = ',', header = 0)
main1=pd.DataFrame(main)
main1.drop("gene_callers_id", axis=1, inplace=True)

main2 = pd.read_csv('ER0043_COG_categories.csv', sep = ',', header = 0) 
main2=pd.DataFrame(main2)
main2.drop("gene_callers_id", axis=1, inplace=True)

merged = pd.merge(main1, main2, left_on = ['source','function','accession'], right_on = ['source','function','accession'], how = 'outer')

#Merge the rest of the samples:

for i in a2:
    main3 = pd.read_csv(''.join(i)+'_COG_categories.csv', sep = ',', header = 0)
    main3 = pd.DataFrame(main3)
    main3.drop('gene_callers_id', axis=1, inplace=True)
    
    merged = pd.merge(merged, main3, left_on = ['source','function','accession'], right_on = ['source','function','accession'], how = 'outer')
    merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(merged)

merged.to_csv('ERIN_NCBI_COG_categories.csv', sep = ',', index = False)


###### Function Frequency ########

# Isolate the frequency of various functions within a sample
  
for i in a1:
    file = pd.read_csv(r''.join(i)+'.ncbi_function.txt', sep='\t', header = None)
    funs = pd.DataFrame(file)
    funs.columns = funs.iloc[0]
    funs.drop(funs.index[0], inplace=True)
    funs.columns = ['gene_callers_id',	'source','accession','function','e_value']
    funs2 = funs[funs.source == 'COG20_FUNCTION']
    freq=pd.DataFrame(funs2['function'].value_counts())
    freq.reset_index(level=0, inplace=True)
    freq.columns = ['function',''.join(i)]
    freq.to_csv(''.join(i)+'_COG_function_frequency.csv', sep = ',', index=False)
    print(i)  
    
    
freq1 = pd.read_csv('ER0003_COG_function_frequency.csv', sep = ',', header = 0) 
freq1=pd.DataFrame(freq1)

freq2 = pd.read_csv('ER0043_COG_function_frequency.csv', sep= ',', header=0)
freq2 = pd.DataFrame(freq2)

freq_merge = pd.merge(freq1, freq2, on='function', how='outer')

for i in a2:
    freq3 = pd.read_csv(''.join(i)+'_COG_function_frequency.csv', sep = ',', header = 0)
    freq3 = pd.DataFrame(freq3)
       
    freq_merge = pd.merge(freq_merge, freq3, on='function', how = 'outer')
    freq_merge.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(freq_merge)

freq_merge.to_csv('ERIN_NCBI_COG_function_frequencies.csv', sep = ',', index = False)


######################################################

# SCG TAXONOMY OUTPUT

######################################################

rootdir = r'D://HPCC/Spring2021_Aim2_Pipeline/Anvio/06_2021_taxonomy_SCG/'
os.chdir(rootdir)


for i in a1:
    file = pd.read_csv(r"".join(i)+'.scg_taxonomy_hits.txt', sep='\t', header = None)
    tax = pd.DataFrame(file)
    tax.columns = tax.iloc[0]
    tax.drop(tax.index[0], inplace=True)
    tax.columns = ['gene_callers_id',	'source','accession','function',''.join(i)]
    tax2 = tax[funs.source == 'COG20_CATEGORY']
    tax2["".join(i)] = pd.to_numeric(funs2["".join(i)])
    tax2.sort_values(by=''.join(i), ascending=True, inplace=True)
    tax2.drop_duplicates(subset ="function", keep = 'first', inplace = True)
#    tax2.to_csv('D://HPCC/Spring2021_Aim2_Pipeline/Anvio/06_2021_function_NCBI_COGS_output/'+"".join(i)+'_COG_categories.csv', sep = ',', index=False)
    print(i)

    file = pd.read_csv(r'ER0003.scg_taxonomy_hits.txt', sep='\t', header = None)
    tax = pd.DataFrame(file)
    tax.columns = tax.iloc[0]
    tax.drop(tax.index[0], inplace=True)
    tax['bitscore'] = pd.to_numeric(tax['bitscore'])
    tax2 = tax[tax.bitscore > 50 ]
    freq=pd.DataFrame(tax2['t_species'].value_counts())
    freq.reset_index(level=0, inplace=True)
    freq.columns = ['species','ER0003']
    freq1=pd.DataFrame(tax2['t_genus'].value_counts())
    freq1.reset_index(level=0, inplace=True)
    freq1.columns = ['genus','ER0003']
    freq2=pd.DataFrame(tax2['t_family'].value_counts())
    freq2.reset_index(level=0, inplace=True)
    freq2.columns = ['family','ER0003']
    freq3=pd.DataFrame(tax2['t_order'].value_counts())
    freq3.reset_index(level=0, inplace=True)
    freq3.columns = ['order','ER0003']
    freq4=pd.DataFrame(tax2['t_class'].value_counts())
    freq4.reset_index(level=0, inplace=True)
    freq4.columns = ['class','ER0003']
    freq5=pd.DataFrame(tax2['t_phylum'].value_counts())
    freq5.reset_index(level=0, inplace=True)
    freq5.columns = ['phylum','ER0003']

#    tax2.to_csv('D://HPCC/Spring2021_Aim2_Pipeline/Anvio/06_2021_function_NCBI_COGS_output/'+"".join(i)+'_COG_categories.csv', sep = ',', index=False)
    print(i)


freq1 = pd.read_csv('ER0003_COG_function_frequency.csv', sep = ',', header = 0) 
freq1=pd.DataFrame(freq1)

freq2 = pd.read_csv('ER0043_COG_function_frequency.csv', sep= ',', header=0)
freq2 = pd.DataFrame(freq2)

freq_merge = pd.merge(freq1, freq2, on='function', how='outer')

for i in a2:
    freq3 = pd.read_csv(''.join(i)+'_COG_function_frequency.csv', sep = ',', header = 0)
    freq3 = pd.DataFrame(freq3)
       
    freq_merge = pd.merge(freq_merge, freq3, on='function', how = 'outer')
    freq_merge.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(freq_merge)

freq_merge.to_csv('ERIN_NCBI_COG_function_frequencies.csv', sep = ',', index = False)



