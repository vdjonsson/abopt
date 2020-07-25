import numpy as np
import pandas as pd

filepath = '../data/'
filename = 'rcsb_pdb_6M0J'

pdb_num = '6M0J'
name_to_find = 'receptor binding domain'
pdb_start = 15 # to change
pdb_end = 208 # to change
offset = 319

with open(filepath+filename+'.fasta', 'r') as file:
    found = False
    for line in file:
        if found:
            sequence = line.strip()
            break
        if name_to_find in line:
            found = True

fasta_pos = np.arange(0, len(sequence)).astype(float)
pdb_pos = np.arange(0, len(sequence)).astype(float)
pdb_pos[np.logical_or(pdb_pos < pdb_start-1, pdb_pos > pdb_end-1)] = np.nan
pdb_pos = pdb_pos + offset

df = pd.DataFrame({'pdb': pdb_num, 'pdb_location': pdb_pos, 'fasta_location': fasta_pos, 'aa_sequence': list(sequence)})

output_filename = 'SARS_CoV_2_RBD'
df.to_csv(filepath+output_filename+'_locations.csv', sep=',', header=True, index=False)
