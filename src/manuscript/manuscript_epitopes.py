import pandas as pd
import seaborn as sb 
import os 

import sys
sys.path.insert(1, '../pipeline/')

import structure as structure
import energy as energy 
import utils as utils 


''' Antibody list '''
abs = utils.get_antibody_data('../../data/meta/antibody_list.txt')

filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'

ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30']
data_dir = '../../data/'
out_dir = '../../output/'

for ab_name in ab_names: 
    
    pdb_name = abs.loc[abs.antibody == ab_name].pdb.values[0]
    repair_dir = out_dir + 'repair/' + pdb_name + '/'
    pdb_dir = data_dir + 'pdb/' + ab_name + '/'

    pdb_file_name = pdb_name + '.pdb'
    labeled_chains = structure.label_chains(pdb_name)

    print (pdb_file_name)
    print (labeled_chains)

    ep = structure.find_epitopes(pdb_dir, pdb_file_name, labeled_chains, 6)    


    print (ep)
    exit()
