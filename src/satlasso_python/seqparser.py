import numpy as np
import pandas as pd

# variables to update:
# filepath (line 35)
# neut_filename (line 36)
# seq_filename (line 37)
# colname (line 44)
# make_sequences args (line 46)

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

def one_hot_encode(aa):
    if aa not in aalist:
        return [0]*len(aalist)
    else:
        encoding = [0]*len(aalist)
        encoding[aalist.index(aa)] = 1
        return encoding

filepath = '../data/'
filename = 'nussenzweig_antibody_data_cleaned_with_alignments'#'single_mut_effects_cleaned'#'kyratsous_neutralization_data'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)

aamatrix = np.empty((0, len(df.sequences[0])*len(aalist)), int)
for sequence in df.sequences:
    row = []
    for aa in sequence:
        row = row + one_hot_encode(aa)
    aamatrix = np.vstack((aamatrix,row))

np.savetxt(filepath+filename+'_sparse_matrix'+'.csv', aamatrix, delimiter=',')
