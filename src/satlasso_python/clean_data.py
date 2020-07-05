import pandas as pd
import numpy as np

def rename_cols(df):
    newcols = []
    for col in df.columns:
        newcol = col
        for ch in [' ', '-']:
            newcol = newcol.replace(ch, '_')
        for ch in ['/','(',')']:
            newcol = newcol.replace(ch, '')
        newcol = newcol.lower()
        newcols.append(newcol)
    return newcols

filepath = '../data/'
filename_neutralization = '41586_2020_2456_MOESM7_ESM'

neut_df = pd.read_csv(filepath+filename_neutralization+'.csv', sep=',', skiprows=2, header=0, nrows=93)

neut_df.columns = rename_cols(neut_df)
neut_df = neut_df.replace('>1000', 1000)
neut_df = neut_df.replace('NT', np.NaN)

filename_sequences = '41586_2020_2456_MOESM6_ESM'

seq_df = pd.read_csv(filepath+filename_sequences+'.csv', sep=',', skiprows=2, header=0)
seq_df = seq_df.loc[:, ~seq_df.columns.str.contains('^Unnamed')]
seq_df.columns = rename_cols(seq_df)

neut_df['igh_vdj_aa'] = seq_df.igh_vdj_aa
neut_df['igl_vj_aa'] = seq_df.igl_vj_aa

output_file = 'nussenzweig_antibody_data'
neut_df.to_csv(filepath+output_file+'_cleaned'+'.csv', header = True, index=False)

filename_bloom = 'single_mut_effects'
df = pd.read_csv(filepath+filename_bloom+'.csv', sep=',', header=0)
df = df[df.mutation.str.contains(r'\*')==False]
df = df.replace('NA', np.NaN)
df = df.dropna()
df.to_csv(filepath+filename_bloom+'_cleaned.csv', sep=',', header=True, index=False)
