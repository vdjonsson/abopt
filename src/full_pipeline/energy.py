import pandas as pd
import os 
import sys 
import foldx as foldx
import structure as su

data_path = '../../data/'
estimator_path = data_path + 'estimator/'
location_path = data_path + 'location/'
pdb_path = data_path + 'pdb/'
out_path = data_path + 'ddg/' 
out_tmp = '../../output/tmp/' 

foldx_path = '/Applications/foldx5MacC11/'
foldx_path = ''

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL','NA']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V' ,'X']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))



def read_dg_file(dgfile):

    df = pd.read_table(dgfile,header=None)
    df = df.rename(columns= {0:'mut_chain', 1:'dg'})
    df = df.drop_duplicates(subset='mut_chain')
    return df


def read_estimator(filename, ab):

    """ Reads filename corresponding to estimator 
    """

    retval = True
    print(filename)
    df = pd.read_csv(filename)
    dfss = df.loc[df['antibody_id'] == ab ]
    if len(dfss) == 0:
        print('No estimator found')
        retval = False
    return dfss 


def calculate_ddg_bind(antibody, pdb, pdb_less, scantype='virus', indir='./', outdir=''):

    less_struct = 'less_ab'
    if scantype == 'ab': less_struct = 'less_virus'

    ddg_path = indir
    prefix='PS_'
    scantype = '_' + scantype
    suffix ='_scanning_output.txt'

    dg1 = read_dg_file(ddg_path + prefix + pdb  + scantype + suffix)
    dg2 = read_dg_file(ddg_path + prefix + pdb_less + scantype + suffix)

    ddg = pd.merge(dg1, dg2, on='mut_chain', suffixes = ['', '_less'])

    # mut = wt, chain, mutation
    ddg['pdb_location'] = ddg.mut_chain.str[4:-1].astype(str)
    ddg['mut'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:]
    ddg['wt'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:4] + ddg.pdb_location
    ddg['chain'] = ddg.mut_chain.str[3]
    ddg['wildtype'] = ddg.mut.str[0] == ddg.mut.str[-1]

    ddg_name = 'ddg_bind'
    ddg[ddg_name] = ddg['dg'] - ddg['dg_less']
    
    return ddg


def constrain_energies(data, cutoffmin =0, cutoffmax=0, topmuts=10):

    data['mutation'] = data.mut.str[-1:]
    data = data.loc[~data.mutation.isin(['e','o'])]
    sorted = data.sort_values(by='ddg_bind').reset_index()
    constrained = sorted.loc[sorted.ddg_bind < cutoffmax ]
    constrained['loc'] = constrained['mut'].str[4:-1]
    constrained = constrained.iloc[0:topmuts+1,:]

    return constrained


def constrain_estimator_features(data, cutoffmin =0, cutoffmax=0, topmuts = 10, filter_name='chain', filter='H'):

    """ Constrain the estimator features 
    """

    indices = (data['wild_type'].astype(bool) == True) & (data['coefficient'].astype(float) >=cutoffmax)
    indices = indices & (data[filter_name] == filter)
    mut = data.loc[indices]
    mut = mut.loc[mut['coefficient'].astype(float).nlargest(n=topmuts).index]
    mut['location'] = mut['location'].astype(str)

    pdb, locations = read_pdb_locations(antibody=data['antibody_id'].values[0])

    merged = mut.merge(locations, how='left', left_on =['chain', 'location'], right_on = ['chain', 'fasta_location'])
    merged['mut'] = merged.wt_pdb + merged.aa_x

    return merged


    
def read_pdb_locations(file_location='', antibody=''):

    """ Read PDB locations for a given file and transform 
    """
    print(file_location)
    if (file_location == ''): 
        file_location = '../../data/location/' + antibody+ '_locations.csv'

    df = pd.read_csv(file_location)
    
    print(df.head())
    df = df.dropna(axis=0)
    df['pdb_location'] = df.pdb_location.astype(str).str[:-2]
    df['fasta_location'] = df.fasta_location.astype(str)
    df['aa_three'] = [aa_1to3_dic[aa] for aa in df.aa.values]
    df['wt_pdb_foldx'] = df.aa + df.chain + df.pdb_location + 'a'
    df['wt_pdb'] = df.aa + df.chain +  df.pdb_location
    pdb = df.pdb.unique()[0]
    df = df.dropna(axis=0)

    return pdb, df


#def calculate_ddgs_from_estimator_coefficients(file_estimator, ab_name, pdb_struct):

#    retval, ab = read_estimator(f,ab_name)
#    if retval == False: 
#        exit()
    #pdb_name, pdb_loc = read_pdb_locations(ab_name)
    #ab_pos = get_mutation_positions(ab)
    
    # pos_scan = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)
    #run_position_scan (ab_name, pdb_struct, pos_scan)
    # now run ddg calculations 


def build_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4):

    mutations = get_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4)
    write_to_mutations_file(ab_name, mutations, 'individual_list.txt')
    run_build_model(ab_name, p + '.pdb', 'individual_list.txt')
    run_build_model(ab_name, p_less + '.pdb', 'individual_list.txt')

    run_multiple_repair_model(ab_name,p, mutations)
    run_multiple_repair_model(ab_name,p_less, mutations)
    return mutations




def write_to_mutations_file(mutations, mutations_filename, out_dir):

    ab_path = pdb_path + ab_name + '/'
    out_path = pdb_path + ab_name + '/'

    mutstr = ''
    for mut in mutations:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutations_filename, 'w')   
    f.write(mutstr)
    f.close()

    f  = open(out_dir+ 'mutations.txt', 'w')
    f.write(mutstr)
    f.close()

    
