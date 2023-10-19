# -*- coding: utf-8 -*-
"""
Nonpareil Output File Merge 

Author: Zoe Hansen
Last Modified: 2022.01.14

This code is designed to take the nonpareil_output.npo files from my external hard drive and 
merge them into a single file for use in R. 

"""
import pandas as pd
import os
from random import randint

#Set root directory
rootdir1=r'D://HPCC/Nonpareil_coverage/'
os.chdir(rootdir1)

# Designate our sample numbers we wish to iterate through
samples = open('ERIN_samples_IDs_clean.txt')
samples1 = samples.read().splitlines()

# Make a table containing our file name, sample ID, and color designations

data=[]
for i in samples1: 
    file = ''.join(i)+'.nonpareil.npo'
    print(file)
    dir_path = rootdir1+file
    print(dir_path)
    data.append((dir_path,i))

npo_df = pd.DataFrame(data, columns = ['File','ER_ID'])
print(npo_df)

# Generating a list of colors 
col_vector = []
n = 270
for i in range(n):
    col_vector.append('#%06X' % randint(0, 0xFFFFFF))


npo_df['color'] = col_vector

print(npo_df)

npo_df.to_csv('Nonpareil_npo_output_files.csv', sep = ',', index = False )





