import numpy as np
import pandas as pd

filepath = '../../data/'
filename = 'NeutSeqData_VH3-53_66'#'NeutSeqData_C002-215_cleaned'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)

colnames = ['VH or VHH', 'VL']#['igh_vdj_aa', 'igl_vj_aa']
identifying_col = 'Name'#'antibody_id'

for colname in colnames:
    with open(filepath+filename+'_'+colname+'.fasta', 'w') as writer:
        for i in range(0,len(df)):
            writer.write('>antibody_'+colname.replace(' ', '_')+'|'+identifying_col+'='+df[identifying_col][i]+'\n')
            writer.write(df[colname][i]+'\n')
