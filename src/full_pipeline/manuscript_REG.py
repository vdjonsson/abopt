import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
import plotutils as pu 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.stats as stats 
from sklearn import preprocessing


filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'

ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30' ,'REG10987', 'REG10933']
pdb_names = ['7bz5', '6xcm', '6xc3', '7c01', '7jmp', '6xe1','6xdg_REG10987', '6xdg_REG10933']
rbd_chains = ['A', 'C', 'C', 'A', 'A', 'E', 'E','E']


abs = ab_names 
pdbs = pdb_names 
ab_names = ['C002-S1', 'C002-S2','C104-S1', 'C110', 'C119', 'C121-S1', 'C121-S2','C135', 'C144']
pdb_names = ['7K8S', '7K8T', '7K8U', '7K8V', '7K8W', '7K8X', '7K8Y', '7K8Z', '7K90']



#abs = abs+ ab_names 
#pdbs = pdbs+ pdb_names 
#ablist = pd.DataFrame()

#ablist['antibody'] = abs 
#ablist['pdb'] = pdbs 


#ablist.to_csv('../../manuscript/antibody_list.txt', index=None)

#exit()


#ab_names = ['B38', 'C105', 'CB6', 'COVA2-39','CV30']
#pdb_names = ['7bz5', '6xcm', '7c01', '7jmp', '6xe1']
#rbd_chains = ['A', 'C', 'A', 'A', 'E']

#ab_names = ab_names[-2:]
#pdb_names = pdb_names[-2:]
#rbd_chains = rbd_chains[-2:]

print (ab_names)
print (pdb_names)

pdb_rbd =dict(zip(pdb_names, rbd_chains))

# input
ab_pdb = dict(zip(ab_names, pdb_names))

pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/' for ab in ab_names]))
repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))

''' Generate all the data: run estimator, energy minimization  '''
'Repair antibody/viral receptor original structure '

print(list(pdb_files.values()))

ap.repair(pdb_dirs = list(pdb_dirs.values()), pdb_list= list(pdb_files.values()), out_dirs=list(repair_dirs.values()))

repaired_wt_pdb = [ pdb +'_Repair.pdb' for pdb in pdb_names]
repaired_wt_dict = dict(zip(pdb_names, repaired_wt_pdb))

' Remove virus and antibody from WT structure and repair these structures '

for pdb_name in pdb_names: 

    ' Remove antibody and virus from original structure and repair '
    labeled_chains = structure.label_chains(pdb_name)
    print(labeled_chains)

    ap.remove(pdb_dirs = [repair_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name]], chains= labeled_chains, chain_type= 'antibody', out_dirs = [remove_dirs[pdb_name]])
    # ap.remove(pdb_dirs = [repair_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name]], chains= labeled_chains, chain_type= 'virus', out_dirs = [remove_dirs[pdb_name]])
    
    ' Repair these structures '
    ap.repair(pdb_dirs = [remove_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name][:-4]+'_less_ab.pdb'], out_dirs =[repair_dirs[pdb_name]])

#    ap.repair(pdb_dirs = [remove_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name][:-4]+'_less_virus.pdb'], out_dirs =[repair_dirs[pdb_name]])

exit()
removed_repaired_pdb = [ pdb +'_Repair_less_ab.pdb' for pdb in pdb_names]
repaired_removed_repaired_pdb = [ pdb +'_Repair_less_ab_Repair.pdb' for pdb in pdb_names]

' Constrain estimator coefficients based on some cutoffs '

for ab_name in ab_names: 

    ' Set up directories '
    p1 = ab_pdb[ab_name] + '_Repair.pdb'
    p2 = ab_pdb[ab_name] + '_Repair_less_ab_Repair.pdb'
 
    ab_list = [p1, p2]
    repair_dir = repair_dirs[ab_pdb[ab_name]]
    scan_dir = scan_dirs[ab_pdb[ab_name]]
    energy_dir = energy_dirs[ab_pdb[ab_name]]
    design_dir = design_dirs[ab_pdb[ab_name]]
    mutate_dir = mutate_dirs[ab_pdb[ab_name]]
    remove_dir = remove_dirs[ab_pdb[ab_name]]
    constrain_dir = constrain_dirs[ab_pdb[ab_name]]    

    lowerbound = 400
    upperbound = 520 

    apdb, df = energy.read_pdb_locations(file_location='../../data/location/SARS_CoV_2_RBD_locations.csv')
    locs = df.loc[(df.pdb_location.astype(int) >= lowerbound) & (df.pdb_location.astype(int) <= upperbound)]                                      
    mutations = locs.aa + pdb_rbd[ab_pdb[ab_name]] + locs.pdb_location +'a'
    scanvalues = mutations.values        

    print(p1)
    print(p1)

    ap.scan (scantype='location', scanvalues = scanvalues, scanmolecule= 'virus', antibody = p1, pdblist = [p1] , pdbdir=repair_dir, outdir=scan_dir)
    ap.scan (scantype='location', scanvalues = scanvalues, scanmolecule= 'virus', antibody = p2, pdblist = [p2], pdbdir=repair_dir, outdir=scan_dir)


