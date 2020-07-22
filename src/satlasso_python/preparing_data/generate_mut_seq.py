import numpy as np
import pandas as pd

# variables to update:
# filepath (line 35)
# neut_filename (line 36)
# seq_filename (line 37)
# colname (line 44)
# make_sequences args (line 46)

def make_sequences(variants, aa_sequence, include_wildtype=False, offset=0):
    variant_sequences = []
    if include_wildtype:
        variant_sequences.append(aa_sequence)
    for variant in variants:
        if variant!='wild_type':
            aa_wild = variant[0]
            index = int(variant[1:len(variant)-1])
            aa_mut = variant[len(variant)-1]
            variant_sequence = list(aa_sequence)
            variant_sequence[index-1+offset] = aa_mut
            variant_sequences.append("".join(variant_sequence))
    return variant_sequences

filepath = '../../data/'
data_filename = 'single_mut_effects_cleaned'#'kyratsous_neutralization_data'
fasta_filename = 'rcsb_pdb_6M0J'

f = open(filepath+fasta_filename+'.fasta', 'r')
f.readline()
aa_sequence = f.readline().strip('\n')
f.close()

colname = 'mutation_RBD' #'variant'
df = pd.read_csv(filepath+data_filename+'.csv', sep=',', header=0)
df['sequences'] = make_sequences(df[colname], aa_sequence, include_wildtype=False, offset=12)

df.to_csv(filepath+data_filename+'_with_sequences.csv', sep=',', header=True, index=False)
