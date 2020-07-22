import numpy as np
import pandas as pd

filepath = '../../data/'
filename = 'nussenzweig_antibody_data_cleaned'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)

author = 'nussenzweig'
colnames = ['igh_vdj_aa', 'igl_vj_aa']
identifying_cols = ['participant_id', 'antibody_id']

for colname in colnames:
    with open(filepath+filename+'_'+colname+'.fasta', 'w') as writer:
        for i in range(0,len(df)):
        writer.write('>antibody_'+colname+'|'+author+'|'+identifying_cols[0]+'='+df[identifying_cols[0]][i]+'|'+identifying_cols[1]+'='+df[identifying_cols[1]][i]+'\n')
            writer.write(df[colname][i]+'\n')

