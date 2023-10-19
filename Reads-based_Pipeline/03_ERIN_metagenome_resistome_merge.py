# -*- coding: utf-8 -*-
"""
@author: Zoe Hansen
Last modified: 2023.09
"""
### This script is designed to merge TSV output files generated within the
### AmrPlusPlus pipeline. The ResistomeAnalyzer and RarefactionAnalyzer output
### contains gene, group, mechanism, and class TSV files which will be merged
### into a single CSV file for all metagenome samples

import pandas as pd
import os
import numpy as np

# Locate deduplicated resistome feature files
rootdir = r'D://HPCC/AmrPlusPlus/Run1_results/SamDedupRunResistome'
os.chdir(rootdir)
    
a=open('Run1_sampleIDs_trimmed.txt') # Sample IDs used to call files in 'for loop'
a1 = a.read().splitlines()

for i in a1:
    gene_file = pd.read_csv(r'D://HPCC/AmrPlusPlus/Run1_results/SamDedupRunResistome/'+"".join(i)+'.gene.tsv', sep='\t', header = None)
    gene_files = pd.DataFrame(gene_file)
    gene_files.columns = gene_files.iloc[0]
    gene_files.drop(gene_files.index[0], inplace=True)
    gene_files.columns = ['Sample','Gene',''+"".join(i), 'Gene_Fraction']
    gene_files.drop("Sample", axis=1, inplace=True)     
#    gene_files.columns = ['Level', 'Value'] # For Rarefaction only 
    gene_files.to_csv('D://HPCC/AmrPlusPlus/Run1_results/SamDedupRunResistome/gene_info_'+"".join(i)+'.csv', sep = ',', index=False)
    print(i)  
  
a.close()


# Merge the resistome gene data (can use this code as a template for the mechanism, group, class data)
# We can use the 'a1' variable to loop through the sample names

# Create a "starter" for the file merge
merged = pd.DataFrame()
c = ['Gene', 'ER0043','ER0073']
a2 = a1[2:]

main = pd.read_csv('gene_info_ER0043.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_27"])
main1=pd.DataFrame(main)
main1.drop('Gene_Fraction', axis=1, inplace=True)
#main1.columns=['Level','ER0043']  #For rarefaction output only

main2 = pd.read_csv('gene_info_ER0138.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_29"])
main2=pd.DataFrame(main2)
main2.drop('Gene_Fraction', axis=1, inplace=True)
#main2.columns=['Level','ER0138']  #For rarefaction output only 

merged = pd.merge(main1, main2, on = 'Gene', how = 'outer')

#Merge the rest of the samples onto the "starter":

for i in a2:
    main3 = pd.read_csv('gene_info_' + ''.join(i) + '.csv', sep = ',', header = 0)
    main3 = pd.DataFrame(main3)
    main3.drop('Gene_Fraction', axis=1, inplace=True)
#    main3.columns=['Level',''+"".join(i)]  #For rarefaction only 
    
    merged = pd.merge(merged, main3, on = 'Gene', how = 'outer')
    merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(merged)

merged.to_csv('Run1_SamDedupRunResistome_Gene_merged.csv', sep = ',', index = False)



#### Repeat this for the Gene fraction Data as well:
merg_frac = pd.DataFrame()
c1 = ['Gene','ER0043','ER0073']

frac = pd.read_csv('gene_info_ER0043.csv', sep=',', header=0)
frac = pd.DataFrame(frac)
frac.drop('ER0043', axis=1, inplace=True)
frac.columns = ['Gene','ER0043_GeneFraction']

frac2 = pd.read_csv('gene_info_ER0073.csv', sep=',', header=0)
frac2 = pd.DataFrame(frac2)
frac2.drop('ER0073', axis=1, inplace=True)
frac2.columns = ['Gene','ER0073_GeneFraction']

merg_frac = pd.merge(frac, frac2, on ='Gene', how='outer')

for i in a2:
    frac3 = pd.read_csv('gene_info_' + ''.join(i) + '.csv', sep = ',', header = 0)
    frac3 = pd.DataFrame(frac3)
    frac3.drop(''+"".join(i), axis=1, inplace=True)
    frac3.columns = ['Gene',''+"".join(i)+'_GeneFraction']
    
    merg_frac = pd.merge(merg_frac, frac3, on = 'Gene', how = 'outer')
    merg_frac.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(merg_frac)

merg_frac.to_csv('Run1_SamDedupResistome_GeneFraction_merged.csv', sep=',', index=False)

### Above code should be repeated for Runs 2-4 (results will be merged at the end of the this code)

###################################################################################

##### Repeat the above code for the Group, Mechanism, and Class level output ###

for i in a1:
    file = pd.read_csv(r'D://HPCC/AmrPlusPlus/Run1_results/SamDedupRunResistome/'+"".join(i)+'.group.tsv',sep = "\t", header = None)
    files = pd.DataFrame(file)
    files.columns = files.iloc[0]
    files.drop(files.index[0], inplace=True)
    files.columns = ['Sample','Group',''+"".join(i)]
    files.drop("Sample", axis=1, inplace=True)    
    files.to_csv('D://HPCC/AmrPlusPlus/Run1_results/SamDedupRunResistome/group_'+"".join(i)+'.csv', sep = ',', index=False)
    print(i)  
  
# Create a "starter" for the file merge
comb = pd.DataFrame()
c = ['Group', 'ER0043','ER0073']
a2 = a1[2:]

main = pd.read_csv('group_ER0043.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_27"])
main1=pd.DataFrame(main)

main2 = pd.read_csv('group_ER0073.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_29"])
main2=pd.DataFrame(main2)

comb = pd.merge(main1, main2, on = 'Group', how = 'outer')
print(comb)
#Merge the rest of the samples:

for i in a2:
    main3 = pd.read_csv('group_' + ''.join(i) + '.csv', sep = ',', header = 0)
    main3 = pd.DataFrame(main3)
    
    comb = pd.merge(comb, main3, on = 'Group', how = 'outer')
    comb.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(comb)

comb.to_csv('Run1_SamDedupResistome_Group_merged.csv', sep = ',', index = False)

###########################################################################################

##### Merge Runs 1-4 into a comprehensive dataframe ####

rootdir=r'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/Second_Analysis/Resistome'
os.chdir(rootdir)

run1 = pd.read_csv('D://HPCC/AmrPlusPlus/Run1_results/RunRarefaction/Run1_SamDedupResistome_Gene_merged.csv', sep=',', header = 0)
run2 = pd.read_csv('D://HPCC/AmrPlusPlus/Run2_results/RunRarefaction/Run2_SamDedupResistome_Gene_merged.csv', sep=',', header = 0)
run3 = pd.read_csv('D://HPCC/AmrPlusPlus/Run3_results/RunRarefaction/Run3_SamDedupResistome_Gene_merged.csv', sep=',', header = 0)
run4 = pd.read_csv('D://HPCC/AmrPlusPlus/Run4_results/RunRarefaction/Run4_SamDedupResistome_Gene_merged.csv', sep=',', header = 0)

merge1 = pd.DataFrame()
merge1 = pd.merge(run1, run2, on = 'Level', how = 'outer')

merge2 = pd.DataFrame()
merge2 =pd.merge(run3, run4, on = 'Level', how = 'outer')

merge_all=pd.DataFrame()
merge_all = pd.merge(merge1, merge2, on = 'Level', how = 'outer')
merge_all.replace(np.nan, 0, inplace=True)

print(merge_all)

merge_all.to_csv('ERIN_metagenomes_SamDedupResistome_Gene.csv', sep=',', index=False)




