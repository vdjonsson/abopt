import pandas as pd 
import seaborn as sb 
import matplotlib.pyplot as plt 





''' This function takes in a sequence and outputs all possible mutations in 
    the string defined in the region vector. 
'''

def mutate_protein(sequence, region): 
    return mutations 


''' Read in list of sequence files and create a dataframe of  sequences 
'''
def read_sequence_files(files, type): 

    sequences = pd.DataFrame()
    for file in files: 
        sequences = pd.read_csv(file)
        sequences = pd.concat([sequences, tmp])
    return sequences



''' The process should go as follows '''

ab_files_dir = '../data/antibody/'
virus_files_dir = '../data/virus/'


''' 
1. Read in antibody and virus bound structure pdb names
2. Read in corresponding fasta file 
3. Read in binding site file (optional) 
4. Mutate every amino acid in every sequence output to mutation input file for FoldX. 
5. Run FoldX command with input file.  Calculate dg on all combinations of mutations on: a) antibody structure unbound b) virus structure unbound c) antibody and virus structure bound 
6. Calculate ddG based on Jonsson thesis derviation. 
```


