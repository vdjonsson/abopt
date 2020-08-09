import pandas as pd
import seaborn as sb 
import os 

data_path = '../../data/'
estimator_path = data_path + 'estimator/'
location_path = data_path + 'location/'
pdb_path = data_path + 'pdb/'

foldx_path = '/Applications/foldx5MacC11/'


aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))


def read_ddg_file(pdb_name, header=None):

    prefix='PS_'
    suffix ='_scanning_output.txt'
    df = pd.read_table('../../data/ddg/' + prefix + pdb_name +suffix,header=header)
    df = df.rename(columns= {0:'mut', 1:'dg'})
    return df

def read_estimator(filename, ab):
    retval = True
    df = pd.read_csv(estimator_path + filename)
    dfss = df.loc[df.antibody_id == ab ]
    if len(dfss) == 0:
        print('No estimator found')
        retval = False
    return retval, dfss 

def calculate_ddg_bind(p1, p2, ddg_name):

    dg1 = read_ddg_file(pdb_name=p1,header=None)
    dg2 = read_ddg_file(pdb_name=p2, header=None)
    ddg = pd.merge(dg1, dg2, on='mut')
    ddg[ddg_name] = ddg['dg_x'] -ddg['dg_y']
    ddg.to_csv('./ddG_' + p1 +'_' + p2 +'.csv')
    return ddg[['mut', ddg_name]]

def construct_positionscan_string (pdb_name, ab_pos, pdb_loc):

    # check to see pdb_name and ab_pos is same 
    check = (ab_pos.antibody_id.unique() == pdb_loc.antibody_id.unique())
    if check == False:
        print('Not matching antibody')
        exit()

    merged = ab_pos.join(pdb_loc, on='fasta_location', how = 'left', lsuffix='_abpos')
    merged = merged.dropna(axis=0, subset=['pdb_wt_foldx'])
    posscan = merged.pdb_wt_foldx.unique()
    pd.Series(posscan).to_csv(data_path + pdb_name + '_position_scan.csv', index=False)
    posscan_str =  ",".join(posscan)
    return posscan_str


# get positions where wildtype== True and coeff>=0, or wildtype==False, and coeff<=0 
def get_mutation_positions(ab):

    # get positions where wildtype== True and coeff>=0 
    # constrain also after reordering 
    mut1 = ab.loc[(ab.wild_type == True) & (ab.coefficient >=1e-10) & (ab.chain == 'H') & (ab.coefficient.nlargest(n=200))]
    mut2 = ab.loc[(ab.wild_type == False) & (ab.coefficient <= -1e-10) & (ab.chain =='H') & (ab.coefficient.nsmallest(n=200))]

    # focus on the heavy chain
    abs = pd.concat([mut1, mut2])
    return abs 
    
def read_pdb_locations(ab_name):
    df = pd.read_csv(location_path + ab_name + '_locations.csv', dtype='object')
    df['pdb_wt_foldx'] = df.aa + df.chain + df.pdb_location + 'a'
    pdb = df.pdb.unique()[0]
    df = df.dropna(axis=0)
    return pdb, df

def run_position_scan (ab_name, pdb_name, pos_scan):

    ab_path = pdb_path + ab_name + '/'
    out_path = data_path + 'ddg/' + ab_name + '/'
    print('Running position scan')
    print(ab_path)
    print (pos_scan)
    command = foldx_path + "foldx --command=PositionScan --pdb-dir=" + ab_path + " --pdb=" + pdb_name
    command = command + " --positions="+ pos_scan +" --out-pdb=false  --output-dir=" + out_path
    print(command)
    os.system(command)

def run_build_model(ab_name, pdb_name, mutations_file): 

    print(pdb_name)
    print(mutations_file)
    command = "foldx --command=BuildModel --pdb="+pdb_name+".pdb "
    command  = command + "--mutant-file=" + mutations_file
    print (command)
    os.system(command)

def run_repair_model(ab_name, pdb_name):
    print(pdb_name)
    command = "foldx --command=RepairPDB --pdb="+pdb_name+".pdb "
    os.system(command)

def run_multiple_repair_model(pdb_name, mutations):
    i = 1
    for mut in muts:
        run_repair_model(pdb_name, pdb_name + '_' + mut)
        i = i+1

def calculate_ddgs_from_estimator_coefficients(file_estimator, ab_name,pdb_struct):
    retval, ab = read_estimator(f,ab_name)
    if retval == False: 
        exit()
    pdb_name, pdb_loc = read_pdb_locations(ab_name)
    ab_pos = get_mutation_positions(ab)
    
    pos_scan = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)
    run_position_scan (ab_name, pdb_struct, pos_scan)


def get_optimal_antibody_structure_mutations(p1, p2, upper_bound_ddg = -0.4):

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
    mutstr = ''
    for mut in muts:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutations_filename, 'w')   
    f.write(mutstr)
    f.close()


# TEA: docstring everything
# TEA:  create a separate main function and file with command line arguments  
# TEA: Estimator file, input as command line argument 
f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'

# TEA: PDB file, input as command line argument 
p1  = '6XCM.pdb'

# TEA: Remove virus from p1, and name xx_less_virus.pdb 
# Repair both structures, output is xx_Repair.pdb and xx_less_virus_Repair.pdb
# run_repair_model(xx.pdb, xx.pdb)
# run_repair_model(xx_less_virus.pdb, xx_less_virus.pdb)

# Calculate ddgs for mutations on these structures, given bound and unbound structures  
# These mutations come from the estimator that Natalie produced. This is used as input to potential mutations


ab_name = 'C105'
p1  = '6XCM_Repair'
f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'

ab_name = 'B38'
p1 = '7bz5_Repair'

ab_name = 'CB6'
p1 = '7c01_Repair'

ab_name = 'CV30'
p1 = '6xe1_Repair'
f = 'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'

ab_name = 'CC12.1'
p1 = '6xc3_CC12.1_Repair'



calculate_ddgs_from_estimator_coefficients(f, ab_name = ab_name, pdb_struct=p1 +'.pdb')
#calculate_ddgs_from_estimator_coefficients(f, ab_name = ab_name, pdb_struct=p1+'_less_virus_Repair.pdb')

# Get the mutations needed to  optimize antibody binding, these are just estimator locations with coefficients < -0.4 

#mutations = get_optimal_antibody_structure_mutations(p1, p2, upper_bound_ddg = -0.4)
#mutations = get_optimal_antibody_structure_mutations(p1, p2, upper_bound_ddg = -0.4)

#write_to_mutations_file(mutations, 'individual_list.txt')

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

#pdb, df = read_pdb_locations('SARS_CoV_2_RBD')
#dfss = df.loc[(df.pdb_location.astype(int) >=400) & (df.pdb_location.astype(int) <= 520)] 
#posscan_str =  ",".join(dfss.pdb_wt_foldx)

# Run position scan on p1_m and p2_m 
#run_position_scan (p1_m, p1_m, posscan_str)
#run_position_scan (p2_m, p2_m, posscan_str)

# Run position scan on p1_wt and p2_wt
p1_wt  = '6XCM_Repair' # repaired wt antibody 
# TEA: remove the antibody from this structure and save to _less_ab.pdb 
p2_wt_nr  = '6XCM_Repair_less_ab'
# Repair 6XCM_Repair_less_ab, output will be 6XCM_Repair_less_ab_Repair.pdb
# run_repair_model(p2_wt_nr, p2_wt_nr)
p2_wt = '6XCM_Repair_less_ab_Repair.pdb'

#run_position_scan (p1_wt, p1_wt, posscan_str)
#run_position_scan (p2_wt, p2_wt, posscan_str)

