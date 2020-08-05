import pandas as pd
import numpy as np

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']
aa_mw = [89, 174, 132, 133, 121, 146, 147, 75, 155, 131, 131, 146, 149, 165, 115,105, 119, 204, 181, 117] # g/mol

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
    
def get_sequence(database, name):
    sequence = database.loc[database.Name == name,:]['VH or VHH'].values[0]+database.loc[database.Name == name,:]['VL'].values[0]
    return sequence
    
def calculate_mw(sequence):
    mw = 0
    for aa in list(sequence):
        mw = mw + aa_mw[aalist.index(aa)]
    return mw
    

filepath = '../../data/'
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

output_file = 'NeutSeqData_C002-215'
neut_df.to_csv(filepath+output_file+'_cleaned'+'.csv', header = True, index=False)

# file = bloom

filename_bloom = 'single_mut_effects'
df = pd.read_csv(filepath+filename_bloom+'.csv', sep=',', header=0)
df = df[df.mutation.str.contains(r'\*')==False]
df = df.replace('NA', np.NaN)
df = df.dropna()
df.to_csv(filepath+filename_bloom+'_cleaned.csv', sep=',', header=True, index=False)

# file = regeneron

db_filename = 'CoV-AbDab_270620'
data_filename = 'NeutSeqData_REGN'

database = pd.read_csv(filepath+db_filename+'.csv', header=0, sep=',')
df = pd.read_csv(filepath+data_filename+'.csv', header=0, sep=',')

for colname in df.columns[1:len(df.columns)]:
    seq = get_sequence(database, colname)
    mw = calculate_mw(seq)
    df[colname] = (df[colname]*mw)*1e6

df.to_csv(filepath+data_filename+'_ngml_units.csv', sep=',', header=True, index=False)
