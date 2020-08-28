import pandas as pd
import numpy as np
import re

filepath = '../../data/'
filename = 'NeutSeqData_C002-215_cleaned'#'NeutSeqData_VH3-53_66'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)
colnames = ['igh_vdj_aa', 'igl_vj_aa']#['VH or VHH', 'VL']
identifying_col = 'antibody_id'#'Name'
clustal_files = {'igh_vdj_aa': 'clustalo-I20200804-202158-0097-73086486-p1m','igl_vj_aa': 'clustalo-I20200804-202404-0723-86608010-p1m'}#{'VH or VHH': 'clustalo-I20200804-205100-0297-30677836-p2m','VL': 'clustalo-I20200804-205215-0756-37061306-p1m'}

for colname in colnames:
    with open(filepath+clustal_files[colname]+'.clustal_num','r') as file:
        antibody_ids = []
        seqs = []
        for line in file:
            if line.startswith('antibody'):
                string_info = list(filter(lambda x: x!='', re.split('  |\t|\n', line)))
                tag = identifying_col+'='
                desc = string_info[0]
                seq = string_info[1].strip(' ')
                antibody_id = desc[desc.find(tag)+len(tag):]
                if antibody_id in antibody_ids:
                    seqs[antibody_ids.index(antibody_id)] = seqs[antibody_ids.index(antibody_id)]+seq
                else:
                    antibody_ids.append(antibody_id)
                    seqs.append(seq)
        aligned_df = pd.DataFrame([antibody_ids, seqs], index=[identifying_col, colname+'_aligned']).T
        df = pd.merge(df, aligned_df)

df['sequence'] = [m+n for m,n in zip(df[colnames[0]+'_aligned'].values,df[colnames[1]+'_aligned'].values)]
df.to_csv(filepath+filename+'_aligned.csv', sep=',', index=False, header=True)
