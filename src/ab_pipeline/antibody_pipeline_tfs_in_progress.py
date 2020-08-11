import os
import pandas as pd
import seaborn as sb
import wget
from biopandas.pdb import PandasPdb

data_path = '../data/'

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR', 'TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T', 'W', 'Y', 'V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))

""" TEA: docstring everything """


def read_ddg_file(pdb_name, header=None):
    """Returns a panda with PS ddG data from pdb"""
    prefix='PS_'
    suffix ='_scanning_output.txt'
    df = pd.read_table(prefix + pdb_name +suffix, header=header, index_col = 0)
    return df


def read_estimator(filename, ab):
    """Returns a panda of data for a specific antibody from a data file"""
    df = pd.read_csv(data_path + filename)
    dfss = df.loc[df.antibody_id == ab]
    return dfss 


def calculate_ddg_bind(p1, p2):
    """given two pdbs with ddgs calcuated, finds ddG of binding, prints to file and returns panda"""
    dg1 = read_ddg_file(pdb_name=p1,header=None)
    dg2 = read_ddg_file(pdb_name=p2, header=None)
    ddg = dg1-dg2
    ddg.to_csv('./ddG_' + p1 +'_' + p2 +'.csv', header=None)
    return ddg


def construct_positionscan_string (pdb_name, ab_pos, pdb_loc):
    """Creates argument string to call PositionScan given pdb"""
    merged = ab_pos.join(pdb_loc, on='fasta_location', how = 'left', lsuffix='_abpos')
    posscan = merged.pdb_wt_foldx.unique()
    pd.Series(posscan).to_csv(data_path + pdb_name + '_position_scan.csv', index=False)
    posscan_str =  ",".join(posscan)
    return posscan_str


def get_mutation_positions(ab):
    """get positions where wildtype== True and coeff>=0, or wildtype==False, and coeff<=0 """
    # get positions where wildtype== True and coeff>=0 
    mut1 = ab.loc[(ab.wild_type == True) & (ab.coefficient >=1e-10)]
    mut2 = ab.loc[(ab.wild_type == False) & (ab.coefficient <= -1e-10)] 
    abs = pd.concat([mut1, mut2])
    return abs 


def read_pdb_locations(ab_name):
    """Reads pdb locations for a given antibody"""
    df = pd.read_csv(data_path + ab_name + '_locations.csv', dtype='object')
    df['pdb_wt_foldx'] = df.aa + df.chain + df.pdb_location + 'a'
    pdb = df.pdb.unique()[0]
    df = df.dropna(axis=0)
    return pdb, df


def run_position_scan(ab_name, pdb_name, pos_scan):
    """runs FoldX PositionScan"""
    print('Running position scan')
    print(pdb_name)
    print (pos_scan)
    command = "foldx --command=PositionScan --pdb="+pdb_name
    command = command + " --positions="+pos_scan +" --out-pdb=false"
    print(command)
    os.system(command)


def run_build_model(ab_name, pdb_name, mutations_file): 
    """runs FoldX BuildModel"""
    print(pdb_name)
    print(mutations_file)
    command = "foldx --command=BuildModel --pdb="+pdb_name+".pdb "
    command  = command + "--mutant-file=" + mutations_file
    print (command)
    os.system(command)


def run_repair_model(ab_name, pdb_name):
    """runs FoldX Repair on pdb"""
    print(pdb_name)
    command = "foldx --command=RepairPDB --pdb="+pdb_name+".pdb "
    os.system(command)


def run_multiple_repair_model(pdb_name, mutations):
    """runs FoldX Repair on all the mutations of pdb"""
    i = 1
    for mut in mutations:
        run_repair_model(pdb_name, pdb_name + '_' + mut)
        i = i+1


def calculate_ddgs_from_estimator_coefficients(file_estimator, ab_name, pdb_struct):
    """read estimator data, run FoldX on interesting ones"""
    ab = read_estimator(f,ab_name)
    pdb_name, pdb_loc = read_pdb_locations(ab_name)
    ab_pos = get_mutation_positions(ab)
    pos_scan = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)
    print (pos_scan)
    run_position_scan (ab_name, pdb_struct, pos_scan)


def get_optimal_antibody_structure_mutations(p1, p2, upper_bound_ddg = -0.4):
    """selects large ddG values"""
    ddg = calculate_ddg_bind(p1, p2)
    ddg_sort = ddg.sort_values(by=1).reset_index()
    ddg_sort = ddg_sort.rename(columns = {0:'mut', 1:'ddg'})
    ddg_less = ddg_sort.loc[ddg_sort.ddg < -0.4 ]
    ddg_less['loc'] = ddg_less['mut'].str[4:-1]
    grouped = ddg_less.groupby(by='loc').count().reset_index().sort_values('ddg', ascending =False).reset_index()
    muts = []
    for location in grouped['loc'].values:
        mut= ddg_less.loc[ddg_less['loc'] == location].iloc[0].mut
        three = mut[0:3]
        tmp = aa_3to1_dic[three] + mut[3:]
        muts.append(tmp)

    return muts


def write_to_mutations_file(mutations, mutations_filename):
    """writes mutations to a file"""
    mutstr = ''
    for mut in mutations:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutations_filename, 'w')   
    f.write(mutstr)
    f.close()



# TEA:  create a separate main function and file with command line arguments  
# TEA: Estimator file, input as command line argument 
f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'

# TEA: PDB file, input as command line argument 
p1  = '6XCM.pdb'


def label_chains(pdb_name):
    fasta_url = 'https://www.rcsb.org/fasta/entry/' + pdb_name
    fasta_file = wget.download(fasta_url)
    out_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == ">":
                broken = line.split("|")
                chain_letters = 1
                chain_name = 2
                chain_id = broken[chain_name]
                out_dict[chain_id] = broken[chain_letters].split(' ')[1].split(',')
    return out_dict


# TEA: Remove virus from p1, and name xx_less_virus.pdb
def rm_virus(pdb_file_name, labeled_chains):
    struc_path = os.path.abspath(data_path + pdb_file_name + ".pdb")
    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    the_pdb = ppdb.df['ATOM']
    lab_chan = [value for key, value in labeled_chains.items() if 'spike' in key.lower()]
    for mol in lab_chan:
        for chain in mol:
            the_pdb = the_pdb.loc[the_pdb['chain_id'] != chain]
    ppdb.df['ATOM'] = the_pdb
    file_path = data_path + pdb_file_name + '_no_virus.pdb'
    ppdb.to_pdb(path=file_path, records=['ATOM'], gz=False, append_newline=True)



    

# Repair both structures, output is xx_Repair.pdb and xx_less_virus_Repair.pdb
# run_repair_model(xx.pdb, xx.pdb)
# run_repair_model(xx_less_virus.pdb, xx_less_virus.pdb)

# Calculate ddgs for mutations on these structures, given bound and unbound structures  
# These mutations come from the estimator that Natalie produced 
p1 = '6XCM_Repair'
p2 = '6XCM_less_virus_Repair'

# calculate_ddgs_from_estimator_coefficients(f, ab_name = 'C105', pdb_struct=p1+'.pdb')
# calculate_ddgs_from_estimator_coefficients(f, ab_name = 'C105', pdb_struct=p2+'.pdb')

# Get the mutations needed to  optimize antibody binding, these are just estimator locations with coefficients < -0.4 
# mutations = get_optimal_antibody_structure_mutations(p1, p2, upper_bound_ddg = -0.4)
# write_to_mutations_file(mutations, 'individual_list.txt')

# Generate mutated PDBs 
# run_build_model('C105' + mut, p1, 'individual_list.txt')

# TEA: rename all the repaired files output from build model from xx_number to xx_mutation

# TEA: remove the antibody from these structures and save to _less_ab.pdb
# TEA: Now repair all mutated structures keep track of these from mutations file 
# run_multiple_repair_model(pdb_name=p1, mutations)

# TEA: For all the repaired structures, (here just as an example of two, but you will get them from the mutations file)
# p1_m  = '6XCM_Repair_TH28D_Repair.pdb'
# p2_m = '6XCM_Repair_TH28D_Repair_less_ab_Repair.pdb'

# run_repair_model(p1_m, p1_m)
# run_repair_model(p2_m, p2_m)

# TEA: find union of ACE2 and AB contact sites given epitopes.csv file 
# TEA: write function to take union of antibody contact site and virus contact site 
# TEA: eg ab loc: 1-4, 8-10; ace2: 2-4, 5  => union of ab+virus contact site: 1-10 

pdb, df = read_pdb_locations('SARS_CoV_2_RBD')
dfss = df.loc[(df.pdb_location.astype(int) >=400) & (df.pdb_location.astype(int) <= 520)] 
posscan_str =  ",".join(dfss.pdb_wt_foldx)

# Run position scan on p1_m and p2_m 
# run_position_scan (p1_m, p1_m, posscan_str)
# run_position_scan (p2_m, p2_m, posscan_str)

# Run position scan on p1_wt and p2_wt
p1_wt  = '6XCM_Repair' # repaired wt antibody 
# TEA: remove the antibody from this structure and save to _less_ab.pdb 
p2_wt_nr  = '6XCM_Repair_less_ab'
# Repair 6XCM_Repair_less_ab, output will be 6XCM_Repair_less_ab_Repair.pdb
# run_repair_model(p2_wt_nr, p2_wt_nr)
p2_wt = '6XCM_Repair_less_ab_Repair.pdb'

run_position_scan (p1_wt, p1_wt, posscan_str)
run_position_scan (p2_wt, p2_wt, posscan_str)

