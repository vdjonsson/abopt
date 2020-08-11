import pandas as pd
import seaborn as sb 
import os 
import sys 

data_path = '../../data/'
estimator_path = data_path + 'estimator/'
location_path = data_path + 'location/'
pdb_path = data_path + 'pdb/'
out_path = data_path + 'ddg/' 

foldx_path = '/Applications/foldx5MacC11/'
foldx_path = ''

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))


def convert_pdb_to_list(ab_name, pdb_name, spec_chain='C' , min_res=400, max_res=520):
    prevnum , type_id = 0,0
    aa,chain,residue_no = 3,4,5
    lis = []
    pdb = pdb_path + ab_name + '/repaired/' + pdb_name +'.pdb' 
    with open(pdb, 'r') as f:
        for line in f:
                broken = line.split(" ")
                broken = list(filter(None, broken))
                if broken[type_id] == "ATOM" and broken[residue_no] != prevnum and broken[residue_no].isnumeric():
                    if not spec_chain or (spec_chain and spec_chain == broken[chain]):
                        for amac in aa_3to1_dic:
                            if amac != broken[aa]:
                                try:
                                    aa_3to1_dic[broken[aa]]
                                except KeyError:
                                    print(broken)
                                    for amac in aa_3to1_dic:
                                        if amac in broken[aa]:
                                            print(amac)
                                            broken[aa] = amac
                                            print(broken)
                        to_list = aa_3to1_dic[broken[aa]] + broken[chain] + broken[residue_no] + 'a'

                        if (int(broken[residue_no]) >= min_res and int(broken[residue_no])<= max_res):
                            lis.append(to_list)
                        
                        prevnum = broken[residue_no]
    return lis

def read_ddg_file(ab_name,pdb_name, header=None, scan_type='virus'):
    ddg_path = out_path + ab_name + '/' 
    prefix='PS_'
    scan_type = '_' + scan_type
    suffix ='_scanning_output.txt'
    df = pd.read_table(ddg_path + prefix + pdb_name +scan_type + suffix,header=header)
    df = df.rename(columns= {0:'mut_chain', 1:'dg'})
    return df

def read_estimator(filename, ab):
    retval = True
    df = pd.read_csv(estimator_path + filename)
    dfss = df.loc[df.antibody_id == ab ]
    if len(dfss) == 0:
        print('No estimator found')
        retval = False
    return retval, dfss 

def calculate_ddg_bind(ab_name,p, p_less, scan_type='virus', mut_name=''):

    ddg_out_path = out_path 
    dg1 = read_ddg_file(ab_name,pdb_name=p,header=None, scan_type=scan_type)
    dg2 = read_ddg_file(ab_name, pdb_name=p_less, header=None, scan_type = scan_type)
    ddg = pd.merge(dg1, dg2, on='mut_chain')
    ddg['mut'] = ddg.mut_chain.str[0:3] + ddg.mut_chain.str[4:]
    ddg_name = 'ddg_' + mut_name
    ddg[ddg_name] = ddg['dg_x'] - ddg['dg_y']
    
    ddg.to_csv(ddg_out_path + 'ddG_' + p + '_' + mut_name  +'.csv')
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

    mut1 = ab.loc[(ab.wild_type == True) & (ab.coefficient.astype(float) >=1e-10) & (ab.chain == 'H')]
    mut1 = mut1.loc[mut1.coefficient.astype(float).nlargest(n=10).index]

    mut2 = ab.loc[(ab.wild_type == False) & (ab.coefficient <= -1e-10) & (ab.chain =='H')]
    mut2 = mut2.loc[mut2.coefficient.astype(float).nsmallest(n=10).index]

    # focus on the heavy chain
    abs = pd.concat([mut1, mut2])

    print (abs)
    return abs 
    
def read_pdb_locations(ab_name):
    df = pd.read_csv(location_path + ab_name + '_locations.csv', dtype='object')
    df['pdb_wt_foldx'] = df.aa + df.chain + df.pdb_location + 'a'
    pdb = df.pdb.unique()[0]
    df = df.dropna(axis=0)
    return pdb, df

def run_position_scan (ab_name, pdb_name, pos_scan):

    ab_path = pdb_path + ab_name + '/repaired/'
    out_dir = out_path + ab_name + '/'

    print('Running position scan')
    print(ab_path)
    print (pos_scan)
    command = foldx_path + "foldx --command=PositionScan --pdb-dir=" + ab_path + " --pdb=" + pdb_name
    command = command + " --positions="+ pos_scan +" --out-pdb=false  --output-dir=" + out_dir
    print(command)
    os.system(command)

def run_build_model(ab_name, pdb_name, mutations_file): 

    ab_path = pdb_path + ab_name + '/'
    out_path = pdb_path + ab_name + '/'
    print(pdb_name)
    command = foldx_path + "foldx --command=BuildModel --pdb-dir=" + ab_path + " --pdb=" + pdb_name
    command = command + " --mutant-file="+ mutations_file + " --output-dir=" + out_path
    print(command)
    os.system(command)

def rename_pdb_files(ab_name, pdb_name, mutations):

    ab_path = pdb_path + ab_name + '/'
    out_path = pdb_path + ab_name + '/'
    print(pdb_name)
    i = 1 
    for mut in mutations:
        pdb =  out_path + pdb_name
        command  = 'mv ' + pdb + '_' + str(i) + '.pdb ' + pdb + '_' + mut + '.pdb'
        print(command)
        os.system(command)
        i = i+1


def run_repair_model(ab_name, pdb_name):

    ab_path = pdb_path + ab_name + '/'
    out_path = pdb_path + ab_name + '/repaired/'
    print(pdb_name)
    command = foldx_path + "foldx --command=RepairPDB --pdb-dir=" + ab_path + " --pdb=" + pdb_name
    command = command +  " --output-dir=" + out_path
    print(command)
    os.system(command)

def run_multiple_repair_model(ab_name, pdb_name, mutations):
    i = 1
    for mut in mutations:
        print(mut)
        run_repair_model(ab_name, pdb_name + '_' + mut +'.pdb')
        i = i+1

def calculate_ddgs_from_estimator_coefficients(file_estimator, ab_name,pdb_struct):
    retval, ab = read_estimator(f,ab_name)
    if retval == False: 
        exit()
    pdb_name, pdb_loc = read_pdb_locations(ab_name)
    ab_pos = get_mutation_positions(ab)
    
    pos_scan = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)
    run_position_scan (ab_name, pdb_struct, pos_scan)

def build_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4):
    mutations = get_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4)
    write_to_mutations_file(ab_name, mutations, 'individual_list.txt')
    run_build_model(ab_name, p + '.pdb', 'individual_list.txt')
    run_build_model(ab_name, p_less + '.pdb', 'individual_list.txt')
    rename_pdb_files(ab_name, p, mutations)
    rename_pdb_files(ab_name, p_less, mutations)
    run_multiple_repair_model(ab_name,p, mutations)
    run_multiple_repair_model(ab_name,p_less, mutations)
    return mutations

def get_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4):

    ddg = calculate_ddg_bind(ab_name,p, p_less)
    print (ddg)
    ddg_sort = ddg.sort_values(by='ddg').reset_index()
    #ddg_sort = ddg_sort.rename(columns = {0:'mut', 1:'ddg'})

    ddg_less = ddg_sort.loc[ddg_sort.ddg < -0.4 ]
    ddg_less['loc'] = ddg_less['mut'].str[4:-1]
    grouped = ddg_less.groupby(by='loc').count().reset_index().sort_values('ddg', ascending =False).reset_index()

    print (ddg_less)
    muts = []
    for location in grouped['loc'].values:
        mut= ddg_less.loc[ddg_less['loc'] == location].iloc[0].mut
        three = mut[0:3]
        tmp = aa_3to1_dic[three] + mut[3:]
        muts.append(tmp)

    return muts


def write_to_mutations_file(ab_name, mutations, mutations_filename):

    ab_path = pdb_path + ab_name + '/'
    out_path = pdb_path + ab_name + '/'

    mutstr = ''
    for mut in mutations:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutations_filename, 'w')   
    f.write(mutstr)
    f.close()

    f  = open(ab_path + 'mutations.txt', 'w')
    f.write(mutstr)
    f.close()

# TEA: docstring everything
# TEA:  create a separate main function and file with command line arguments  
# TEA: Estimator file, input as command line argument 
#f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'

# TEA: PDB file, input as command line argument 
#p1  = '6XCM.pdb'

# TEA: Remove virus from p1, and name xx_less_virus.pdb 
# Repair both structures, output is xx_Repair.pdb and xx_less_virus_Repair.pdb
# run_repair_model(xx.pdb, xx.pdb)
# run_repair_model(xx_less_virus.pdb, xx_less_virus.pdb)

# Calculate ddgs for mutations on these structures, given bound and unbound structures  
# These mutations come from the estimator that Natalie produced. This is used as input to potential mutations

#ab_name = 'C105'
#p1  = '6XCM_Repair'
#f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'

#ab_name = 'B38'
#p1 = '7bz5_Repair'

ab_name = 'CB6'
p = '7c01_Repair'

#ab_name = 'CV30'
#p1 = '6xe1_Repair'
#f = 'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'

# doesnt' work i don't know why 
#ab_name = 'CC12.1'
#p1 = '6xc3_CC12.1_Repair'


#calculate_ddgs_from_estimator_coefficients(f, ab_name = ab_name, pdb_struct=p1 +'.pdb')
# Get the mutations needed to  optimize antibody binding, these are just estimator locations with coefficients < -0.4 

less_str = '_less_virus_Repair'

#build_optimal_antibody_structure_mutations(ab_name,p, p+less_str, upper_bound_ddg = -0.4)

#run_multiple_repair_model(pdb_name=p, mutations)


# Generate mutated PDBs 
# run_build_model('C105' + mut, p1, 'individual_list.txt')

# TEA: rename all the repaired files output from build model from xx_number to xx_mutation 
# TEA: remove the antibody from these structures and save to _less_ab.pdb
# TEA: Now repair all mutated structures keep track of these from mutations file 

#run_multiple_repair_model(pdb_name=p1, mutations)
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
#p1_wt  = '6XCM_Repair' # repaired wt antibody 
# TEA: remove the antibody from this structure and save to _less_ab.pdb 
#p2_wt_nr  = '6XCM_Repair_less_ab'
# Repair 6XCM_Repair_less_ab, output will be 6XCM_Repair_less_ab_Repair.pdb
# run_repair_model(p2_wt_nr, p2_wt_nr)
#p2_wt = '6XCM_Repair_less_ab_Repair.pdb'

#run_position_scan (p1_wt, p1_wt, posscan_str)
#run_position_scan (p2_wt, p2_wt, posscan_str)

if __name__ == "__main__":
    args = sys.argv
    ab_name = args[1]
    pdb_name = args[2]
    chain_name = args[3]

    less_str = '_less_virus_Repair'

    f = 'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
    #calculate_ddgs_from_estimator_coefficients(f, ab_name = ab_name, pdb_struct=pdb_name +'.pdb')
    #build_optimal_antibody_structure_mutations(ab_name,pdb_name, pdb_name+less_str, upper_bound_ddg = -0.4)
    
    pos_list = convert_pdb_to_list(ab_name, pdb_name, spec_chain=chain_name)
    pos_scan =  ",".join(pos_list)
    run_position_scan (ab_name, pdb_name + '.pdb', pos_scan)
    
    #run_pos_scan(ab_name, pdb_name, convert_pdb_to_list(ab_name, pdb_name, spec_chain=chain_name))
    
