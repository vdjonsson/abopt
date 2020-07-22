import pandas as pd
import numpy as np
import re

filepath = '../../data/'
filename = 'nussenzweig_antibody_data_cleaned'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)
colnames = ['igh_vdj_aa', 'igl_vj_aa']
clustal_files = {'igh_vdj_aa': 'clustalo-I20200702-233207-0599-10860100-p2m','igl_vj_aa': 'clustalo-I20200702-232723-0835-60033944-p2m'}

for colname in colnames:
    with open(filepath+clustal_files[colname]+'.clustal_num','r') as file:
        antibody_ids = []
        seqs = []
        for line in file:
            if line.startswith('antibody'):
                string_info = list(filter(lambda x: x!='', re.split('  |\t|\n', line)))
                tag = 'antibody_id='
                desc = string_info[0]
                seq = string_info[1].strip(' ')
                antibody_id = desc[desc.find(tag)+len(tag):]
                if antibody_id in antibody_ids:
                    seqs[antibody_ids.index(antibody_id)] = seqs[antibody_ids.index(antibody_id)]+seq
                else:
                    antibody_ids.append(antibody_id)
                    seqs.append(seq)
        aligned_df = pd.DataFrame([antibody_ids, seqs], index=['antibody_id', colname+'_aligned']).T
        aligned_df = aligned_df.sort_values(by=['antibody_id'])
        df[colname+'_aligned']=aligned_df[colname+'_aligned'].values

df['sequences'] = [m+n for m,n in zip(df[colnames[0]+'_aligned'].values,df[colnames[1]+'_aligned'].values)]
df.to_csv(filepath+filename+'_with_alignments.csv', sep=',', index=False, header=True)
