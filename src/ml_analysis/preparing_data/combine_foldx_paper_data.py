import pandas as pd
import numpy as np
import csv

filepath = '../../data/'
foldx_filename = '6vw1'
data_filename = 'single_mut_effects_cleaned'

df = pd.read_csv(filepath+data_filename+'.csv', sep=',', header=0)
df['ddG'] = [np.NaN]*len(df)

with open(filepath+foldx_filename+'.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[0][1] == 'E':
            mutation = row[0].split(' ')[0][:1] + row[0].split(' ')[0][2:]
            ddG = row[0].split(' ')[1]
            df.ddG[df.mutation == mutation] = ddG

for i in range(0,len(df.mutation)):
    if df.mutation[i][0] == df.mutation[i][len(df.mutation[i])-1]:
        df.ddG[i] = 0

cols_to_include = ['mutation', 'site_SARS2', 'mutant', 'bind_avg', 'ddG']
df[cols_to_include].to_csv(filepath+'combined_'+data_filename+'_foldx_'+foldx_filename+'.csv', sep=',', header=True, index=False)
