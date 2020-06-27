import pandas as pd 
import seaborn as sb 
import matplotlib.pyplot as plt 


''' This function takes in a sequence and outputs all possible mutations in 
    the string defined in the region vector. 
'''
aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL','STOP']
aa_single = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'R','V','*']
aa_dict = dict(zip(aa_three, aa_single))


def mutate_protein(sequence_name): 
    print(aa_dict)

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

sequence_name ='aPDB'
mutate_protein(sequence_name)

''' 
1. Read in antibody and virus bound structure pdb names
2. Read in corresponding fasta file 
3. Read in binding site file (optional) 
4. Mutate every amino acid in every sequence output to mutation input file for FoldX. 
5. Run FoldX command with input file.  Calculate dg on all combinations of mutations on: a) antibody structure unbound b) virus structure unbound c) antibody and virus structure bound 
6. Calculate ddG based on Jonsson thesis derviation. 
'''
