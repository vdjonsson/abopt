import pandas as pd
import seaborn as sb 
import os 
import sys 
import antibody_pipeline as ap
import energy     


def run_virus_scan(ab):

    ''' Read ab configuration file '''
    abdf = pd.read_table('../../manuscript/antibody_list.txt', sep=',') 
    
    ab_names = abdf.antibody.astype(str).values
    pdb_names = abdf.pdb.astype(str).values
    rbd_chains = abdf.rbdchain.astype(str).values

    print(ab_names) 
    pdb_rbd =dict(zip(pdb_names, rbd_chains))
    ab_pdb = dict(zip(ab_names, pdb_names))


    pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
    pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/unrepaired/' for ab in ab_names]))
    repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
    remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
    constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
    scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
    energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
    design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
    mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))


    ab_name = ab
    pdb_name = ab_pdb[ab]
    p1 = pdb_name + '_Repair.pdb'
    p2 = pdb_name + '_Repair_less_ab_Repair.pdb'
    
    ab_list = [p1, p2]

    print(ab_list) 

    repair_dir = repair_dirs[ab_pdb[ab_name]]
    scan_dir = scan_dirs[ab_pdb[ab_name]]
    
    ''' Run mutational scanning on viral receptor unbound and bound to mutated antibody '''
    lowerbound = 400
    upperbound = 400

    print('about to read energy location')
    apdb, df = energy.read_pdb_locations(file_location='../../data/location/SARS_CoV_2_RBD_locations.csv')

    print(apdb)
    print(df)

    locs = df.loc[(df.pdb_location.astype(int) >= lowerbound) & (df.pdb_location.astype(int) <= upperbound)]                 
    print(locs.pdb_location.astype(str) + 'a')
    print(pdb_rbd[ab_pdb[ab_name]])
    print(locs.aa.dtype)
    mutations = locs.aa + pdb_rbd[ab_pdb[ab_name]] + locs.pdb_location.astype(str) +'a'
    print(mutations)
    scanvalues = mutations.values    
    print(scanvalues)
    print('about to scan')
    ap.scan (scantype='location', scanvalues = scanvalues, scanmolecule= 'virus', antibody = ab_name, pdblist = ab_list , pdbdir=repair_dir, outdir=scan_dir)


if __name__ == "__main__":

    ab = ''
    print(f"Arguments count: {len(sys.argv)}")
    for i, arg in enumerate(sys.argv):
        ab = sys.argv[1]
        print (ab) 
        
    run_virus_scan(ab)

