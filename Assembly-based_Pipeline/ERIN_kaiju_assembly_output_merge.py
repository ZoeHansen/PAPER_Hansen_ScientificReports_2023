# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 13:37:30 2021

@author: Zoe Hansen
Last Modified: 2021.08.21
"""

### Similar to other 'merge' scripts, this code will take the output from the 
### 'Kaiju_assembly_annotation_072021' folder and concatenate all sample information together

import pandas as pd
import os
import numpy as np

rootdir = r'D://HPCC/Spring2021_Aim2_Pipeline/Kaiju_assembly_annotation_072021/'
os.chdir(rootdir)
    
a=open('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ERIN_samples_IDs_clean.txt')
a1 = a.read().splitlines()

for i in a1:
    file = pd.read_csv(r"".join(i)+'.kaiju_genecalls_summary_all.tsv', sep='\t', header = None)
    mods = pd.DataFrame(file)
    mods.columns = mods.iloc[0]
    mods.drop(mods.index[0], inplace=True)
    mods[['Domain','Phylum','Class','Order','Family','Genus','Species']] = mods.taxon_name.str.split(sep=';', expand=True)
    mods1 = mods[['module_name','module_category','module_subcategory','module_completeness']]
    mods1["module_completeness"] = pd.to_numeric(mods1["module_completeness"])
#    mods.drop(['unique_id','contig_name','kegg_module','module_definition','kofam_hits_in_module','gene_caller_ids_in_module'], axis=1, inplace=True)
    mods2 = mods1[mods1.module_completeness > 0.70]
    mods2.columns = ['module_name','module_category','module_subcategory',''+"".join(i)]
    mods2.to_csv('D://HPCC/Spring2021_Aim2_Pipeline/Anvio/072021_KEGG_metabolism_annotations/'+"".join(i)+'KEGG_module_70cutoff.csv', sep = ',', index=False)
    print(i)  
  
a.close()


    file = pd.read_csv(r'ER0003.kaiju_genecalls_summary_all.tsv', sep='\t', header = None)
    mods = pd.DataFrame(file)
    mods.columns = mods.iloc[0]
    mods.drop(mods.index[0], inplace=True)
    mods[['Domain','Phylum','Class','Order','Family','Genus','Species','random']] = mods.taxon_name.str.split(pat=';', expand=True)
    mods.drop(['file','taxon_id','taxon_name','random'], axis=1, inplace=True)
   
#### Figure out how to "groupby" different levels and add percents and reads (separate dataframes)

    
    freq=pd.DataFrame(mods['Species'].value_counts())
    freq.reset_index(level=0, inplace=True)
    freq.columns = ['species','ER0003']
    freq1=pd.DataFrame(mods['Genus'].value_counts())
    freq1.reset_index(level=0, inplace=True)
    freq1.columns = ['genus','ER0003']
    freq2=pd.DataFrame(mods['Family'].value_counts())
    freq2.reset_index(level=0, inplace=True)
    freq2.columns = ['family','ER0003']
    freq3=pd.DataFrame(mods['Order'].value_counts())
    freq3.reset_index(level=0, inplace=True)
    freq3.columns = ['order','ER0003']
    freq4=pd.DataFrame(mods['Class'].value_counts())
    freq4.reset_index(level=0, inplace=True)
    freq4.columns = ['class','ER0003']
    freq5=pd.DataFrame(mods['Phylum'].value_counts())
    freq5.reset_index(level=0, inplace=True)
    freq5.columns = ['phylum','ER0003']


# Merge the resistome gene data (can use this code as a template for the mechanism, group, class data)
# We can use the 'a1' variable to loop through the sample names

# Create a "starter" for the file merge
merged = pd.DataFrame()
c = ['module_name','module_category','module_subcategory', 'ER0003','ER0043']
a2 = a1[2:]

main = pd.read_csv('ER0003KEGG_module_70cutoff.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_27"])
main1=pd.DataFrame(main)
main1.sort_values(by='ER0003', ascending=False, inplace=True)
main1.drop_duplicates(subset ="module_name", keep = 'first', inplace = True)

main2 = pd.read_csv('ER0043KEGG_module_70cutoff.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_29"])
main2=pd.DataFrame(main2)
main2.sort_values(by='ER0043', ascending=False, inplace=True)
main2.drop_duplicates(subset='module_name', keep='first', inplace=True)

merged = pd.merge(main1, main2, left_on = ['module_name','module_category','module_subcategory'], right_on = ['module_name','module_category','module_subcategory'], how = 'outer')

merged.to_csv('D://HPCC/Spring2021_Aim2_Pipeline/Anvio/072021_KEGG_metabolism_annotations/02_test.csv', sep=',', index=False)