import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     

filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'

ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30']
pdb_names = ['7bz5', '6xcm', '6xc3', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'C', 'A', 'A', 'E']

ab_names = ['B38', 'C105', 'CB6', 'COVA2-39','CV30']
pdb_names = ['7bz5', '6xcm', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'A', 'A', 'E']

ab_names = ab_names[3:]
pdb_names = pdb_names[3:]
rbd_chains = rbd_chains[3:]

print (ab_names)
print (pdb_names)

pdb_rbd =dict(zip(pdb_names, rbd_chains))

# input

ab_pdb = dict(zip(ab_names, pdb_names))

print(ab_pdb)
pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/' for ab in ab_names]))
repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))


v = pd.Series()
for ab_name in ab_names: 

    pdb_name = ab_pdb[ab_name]
    repair_dir = repair_dirs[pdb_name]

    pdb_file_name = pdb_name + '.pdb'
    labeled_chains = structure.label_chains(pdb_name)
    ep = structure.find_epitopes(pdb_dirs[pdb_name], pdb_file_name, labeled_chains, 6)    

    epitope = pd.read_csv('../../data/location/' + pdb_name + '_epitopes.csv')

    design = pd.read_csv('../../output/design/' + ab_name + '/' + ab_name + '_design.csv')

    mutations = design.mut_x
    locations = design.pdb_location_x.unique()

    chain = pdb_rbd[pdb_name]
    virus_locations = epitope.loc[epitope.number_antibody.isin(locations) & (epitope.chain_virus == chain)].number_virus.values

    remove_dir = remove_dirs[pdb_name]

    for mut in mutations: 
        print ('mutation ' + mut)
        mut_file = pdb_name + '_Repair_' + mut + '_Repair.pdb'
        mut_file_less_ab = pdb_name + '_Repair_' + mut + '_Repair_less_ab.pdb'

        ap.repair(pdb_dirs = [remove_dir], pdb_list=[mut_file_less_ab], out_dirs=[repair_dir])

        #exit()
        #print(mut_file)
        #print(mut_file_less_ab)
        #exit()
        #for v in virus_locations:
        #    print( 'virus location:' +str(v))


    



    
