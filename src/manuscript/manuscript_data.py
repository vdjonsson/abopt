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


ab_names = ab_names[1:]
pdb_names = pdb_names[1:]
rbd_chains = rbd_chains[1:]

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

''' Generate all the data: run estimator, energy minimization  '''



'Repair antibody/viral receptor original structure '

#ap.repair(pdb_dirs = list(pdb_dirs.values()), pdb_list= list(pdb_files.values()), out_dirs=list(repair_dirs.values()))

repaired_wt_pdb = [ pdb +'_Repair.pdb' for pdb in pdb_names]
repaired_wt_dict = dict(zip(pdb_names, repaired_wt_pdb))

' Remove virus and antibody from WT structure and repair these structures '

'''
for pdb_name in pdb_names: 

    ' Remove antibody and virus from original structure and repair '
    labeled_chains = structure.label_chains(pdb_name)
    ap.remove(pdb_dirs = [repair_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name]], chains= labeled_chains, chain_type= 'antibody', out_dirs = [remove_dirs[pdb_name]])
    ap.remove(pdb_dirs = [repair_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name]], chains= labeled_chains, chain_type= 'virus', out_dirs = [remove_dirs[pdb_name]])

    ' Repair these structures '
    ap.repair(pdb_dirs = [remove_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name][:-4]+'_less_ab.pdb'], out_dirs =[repair_dirs[pdb_name]])
    ap.repair(pdb_dirs = [remove_dirs[pdb_name]], pdb_list = [repaired_wt_dict[pdb_name][:-4]+'_less_virus.pdb'], out_dirs =[repair_dirs[pdb_name]])

'''

removed_repaired_pdb = [ pdb +'_Repair_less_virus.pdb' for pdb in pdb_names]
repaired_removed_repaired_pdb = [ pdb +'_Repair_less_virus_Repair.pdb' for pdb in pdb_names]

' Constrain estimator coefficients based on some cutoffs '

for ab_name in ab_names: 
    print(ab_name)

    constrain_dir = constrain_dirs[ab_pdb[ab_name]]    
    
    print(constrain_dir)
    #if  'CC121' in ab_name:
    #    ab_name ='CC12.1'
        
    # ap.constrain(constraintype ='estimator', constrainfile=filename, antibody=ab_name, cutoff = [-0.4, 0.1], top=10,out_dir = constrain_dir)

    ' Scan the antibody at locations found by estimator '

    file_estimator = '../../output/constrain/' + ab_name + '/' + ab_name + '_estimator.csv'
    estimator = pd.read_csv(file_estimator).wt_pdb_foldx

    print(estimator)
    
    p1 = ab_pdb[ab_name] + '_Repair.pdb'
    p2 = ab_pdb[ab_name] + '_Repair_less_virus_Repair.pdb'
 
    ab_list = [p1, p2]
    repair_dir = repair_dirs[ab_pdb[ab_name]]
    scan_dir = scan_dirs[ab_pdb[ab_name]]
    energy_dir = energy_dirs[ab_pdb[ab_name]]
    design_dir = design_dirs[ab_pdb[ab_name]]
    mutate_dir = mutate_dirs[ab_pdb[ab_name]]
    remove_dir = remove_dirs[ab_pdb[ab_name]]

    '''
    for p in ab_list:
        ap.scan (scantype='location', scanvalues = estimator, scanmolecule= 'ab', antibody = ab_name, pdblist = [p], pdbdir=repair_dir, outdir=scan_dir)
    '''

    # ap.energy (antibody=ab_name, pdb=p1[:-4], pdb_less=p2[:-4], scantype = 'ab', energy_type='ddgbind', indir = scan_dir, outdir=energy_dir)

    ' Constrain mutation locations where ddG < 0, wild type == True, based on some cutoffs '
    filename= energy_dir + 'ddgbind_' + ab_name + '_ab_scanning.txt'

    #ap.constrain (constraintype ='energy', constrainfile=filename, antibody=ab_name, cutoff = [-1e3, 0.4], top=10, out_dir=constrain_dir)
    #ap.design (designtype='antibody design', designval = ab_name ,file_estimator= constrain_dir + ab_name +'_estimator.csv', file_energies= constrain_dir + ab_name +'_energies.csv', out_dir=design_dir)

    ' Mutate antibody and save the new PDB structure '
    mutations = pd.read_csv(design_dir + ab_name + '_design.csv').mut_x

    # ap.mutate(p1, mutations, repair_dir, mutate_dir)
    mutated_pdbs = [ p1[:-4] + '_' + mut + '.pdb' for mut in mutations]

    for mutated_pdb in mutated_pdbs:

        print (mutated_pdb)
        ' Repair the mutated antibody and save the new PDB structure '
        ap.repair(pdb_dirs = [mutate_dir], pdb_list=[mutated_pdb], out_dirs=[repair_dir])

        repaired_mutated_pdb = mutated_pdb[:-4] + '_Repair.pdb'

        ' Remove antibody from structure and repair the unbound viral receptor, save to new PDB structure ' 
        print(mutated_pdb[0:4])
        labeled_chains = structure.label_chains(mutated_pdb[0:4])

        ap.remove(pdb_dirs = [repair_dir], pdb_list = [repaired_mutated_pdb] , chains= labeled_chains, chain_type= 'antibody', out_dirs = [remove_dir])
        ap.repair(pdb_dirs = [remove_dir], pdb_list=[repaired_mutated_pdb[:-4] + '_less_ab.pdb'], out_dirs=[repair_dir])                                       
        
    exit()
    
    ''' Run mutational scanning on viral receptor unbound and bound to mutated antibody '''

    ''' Constrain the locations for viral scanning '''
    
    # specify chain for the structure as well, TEA find a way to do this with biopandas easily 
    '''
    lowerbound = 400
    upperbound = 400 #520 

    apdb, df = energy.read_pdb_locations(file_location='../../data/location/SARS_CoV_2_RBD_locations.csv')
    locs = df.loc[(df.pdb_location.astype(int) >= lowerbound) & (df.pdb_location.astype(int) <= upperbound)]                                      
    repaired_removed_repaired_mutated_pdbs = [mut[:-4] +  '_Repair.pdb' for mut in removed_repaired_mutated_pdbs]
    
    repaired_wt_pdb = [pdb_name + '_Repair.pdb' ]
    repaired_removed_wt_pdb = [ pdb_name + '_Repair_less_ab_Repair.pdb' ]
    
    pdblist = [repaired_mutated_pdbs[0]] + [repaired_removed_repaired_mutated_pdbs[0]] + repaired_wt_pdb + repaired_removed_wt_pdb
    
    mutations = locs.aa + pdb_rbd[pdb_name] + locs.pdb_location +'a'
    scanvalues = mutations.values    
    
    print(pdblist)
    print(mutations)
    
    ap.scan (scantype='location', scanvalues = scanvalues, scanmolecule= 'virus', antibody ='C105', pdblist = pdblist , pdbdir=repair_dir, outdir=scan_dir)
    '''

    ''' With the epitopes file find the epitope of the mutated antibody on RBD , eg: for H27S mutation on antibody, (chain H, loc 27) find RBD proximal locations  ''' 

    # epitopes = structure.find_epitopes(pdb_file_name, labeled_chains)
    
    ''' Calculate the dddG antibody/ viral mutations  '''
