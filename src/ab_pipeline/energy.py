import pandas as pd
import os 
import sys 


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



def read_dg_file(ab_name,pdb_name, header=None, scantype='virus'):
    ddg_path = out_path + ab_name + '/' 
    prefix='PS_'
    scantype = '_' + scantype
    suffix ='_scanning_output.txt'
    df = pd.read_table(ddg_path + prefix + pdb_name +scantype + suffix,header=header)
    df = df.rename(columns= {0:'mut_chain', 1:'dg'})
    df = df.drop_duplicates(subset='mut_chain')
    return df

def read_estimator(filename, ab):
    retval = True
    df = pd.read_csv(estimator_path + filename)
    dfss = df.loc[df.antibody_id == ab ]
    if len(dfss) == 0:
        print('No estimator found')
        retval = False
    return retval, dfss 

def calculate_ddg_bind(ab_name,p, p_less, scantype='virus', mut_name=''):

    less_struct = 'ab'
    if scantype == 'ab': less_struct = 'virus'
    ddg_out_path = out_path 
    dg1 = read_ddg_file(ab_name,pdb_name=p,header=None, scantype=scantype)
    dg2 = read_ddg_file(ab_name, pdb_name=p_less, header=None, scantype = scantype)

    ddg = pd.merge(dg1, dg2, on='mut_chain', suffixes = ['', '_less'])

    print(ddg.head())
    
    # mut = wt, chain, mutation
    ddg['pdb_location'] = ddg.mut_chain.str[4:-1].astype(str)
    ddg['mut'] = convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:]
    ddg['wt'] = convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:4] + ddg.pdb_location
    ddg['chain'] = ddg.mut_chain.str[3]
    ddg['wildtype'] = ddg.mut.str[0] == ddg.mut.str[-1]

    ddg_name = 'ddg_bind'
    ddg[ddg_name] = ddg['dg'] - ddg['dg_less']
    
    ddg.to_csv(ddg_out_path + 'ddg_bind_' + ab_name + '_' + p + '_' + scantype + '_scanning' + mut_name  +'.csv')
    return ddg


# get positions where wildtype== True and coeff>=0, or wildtype==False, and coeff<=0 and return with pdb locations
def get_antibody_mutation_positions(data, cutoffmin =0, cutoffmax=0, topmuts = 10, filter_name='chain', filter='H'):

    # get positions where wildtype== True and coeff > cutoffmax 
    mut1 = data.loc[(data.wild_type == True) & (data.coefficient.astype(float) >=cutoffmax) & (data[filter_name] == filter)]
    mut1 = mut1.loc[mut1.coefficient.astype(float).nlargest(n=topmuts).index]

    # get positions where wildtype==False and coeff < -cutoffmin
    mut2 = data.loc[(data.wild_type == False) & (data.coefficient <= cutoffmin) & (data[filter_name] == filter)]
    mut2 = mut2.loc[mut2.coefficient.astype(float).nsmallest(n=topmuts).index]

    estimator = pd.concat([mut1, mut2])
    
    # convert estimator fasta locations to pdb locations
    pdb, locations = read_pdb_locations(data.Name.values[0])

    # merge estimator with pdb locations 
    merged = estimator.merge(locations, how='left', left_on =['chain', 'fasta_location'], right_on = ['chain', 'fasta_location'])
    print(merged.head())
    merged['mut'] = merged.wt_pdb + merged.aa_x

    return merged

    
def read_pdb_locations(ab_name):
    df = pd.read_csv(location_path + ab_name + '_locations.csv')
    df['pdb_location'] = df.pdb_location.astype(str).str[:-2]
    df['fasta_location'] = df.fasta_location.astype(str)
    df['aa_three'] = [aa_1to3_dic[aa] for aa in df.aa.values]
    df['wt_pdb_foldx'] = df.aa_three + df.chain + df.pdb_location + 'a'
    df['wt_pdb'] = df.aa + df.chain +  df.pdb_location
    pdb = df.pdb.unique()[0]
    df = df.dropna(axis=0)
    return pdb, df


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

def get_optimal_antibody_structure_mutations(ab_name, pdb, pdb_less, upper_bound_ddg = -0.4):

    # read file instead here 
    ddg = calculate_ddg_bind(ab_name,pdb, pdb_less)
    ddg_sort = ddg.sort_values(by='ddg_bind').reset_index()

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

    
