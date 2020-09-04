import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
#import structure 
import energy     

filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'

ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30']
pdb_names = ['7bz5', '6xcm', '6xc3', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'C', 'A', 'A', 'E']

ab_names = ['B38', 'C105', 'CB6', 'COVA2-39','CV30']
pdb_names = ['7bz5', '6xcm', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'A', 'A', 'E']


ab_names = [ab_names[0]]
pdb_names = [pdb_names[0]]
rbd_chains = [rbd_chains[0]]

print (ab_names)
print (pdb_names)

pdb_rbd =dict(zip(pdb_names, rbd_chains))

# input

ab_pdb = dict(zip(ab_names, pdb_names))

print(ab_pdb)
pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/unrepaired/' for ab in ab_names]))
repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))



ab_name = 'B38'
p1 = '7bz5_Repair.pdb'
p2 = '7bz5_Repair_less_ab_Repair.pdb'

ab_list = [p1, p2]
repair_dir = repair_dirs[ab_pdb[ab_name]]
scan_dir = scan_dirs[ab_pdb[ab_name]]

''' Run mutational scanning on viral receptor unbound and bound to mutated antibody '''

lowerbound = 400
upperbound = 520 

print('about to read energy location')
apdb, df = energy.read_pdb_locations(file_location='../../data/location/SARS_CoV_2_RBD_locations.csv')
locs = df.loc[(df.pdb_location.astype(int) >= lowerbound) & (df.pdb_location.astype(int) <= upperbound)]                                      
mutations = locs.aa + pdb_rbd[ab_pdb[ab_name]] + locs.pdb_location +'a'
scanvalues = mutations.values    
    
print('about to scan')
ap.scan (scantype='location', scanvalues = scanvalues, scanmolecule= 'virus', antibody = ab_name, pdblist = pdblist , pdbdir=repair_dir, outdir=scan_dir)
