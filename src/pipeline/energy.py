import pandas as pd
import os 
import sys 
import structure as su
import matplotlib.pyplot as plt 
import seaborn as sb 
import utils as utils
import argparse

'''
### energy        
`abopt energy` generates a matrix of folding energies based on running energy minimization for mutational scanning. The arguments for the energy command are:

> `input: binding str to calculate binding ddG `

> `input: coupling str to calculate coupling dddGs `

> `input: locations list of locations to constrain energy calculations, comma delimited, eg: 250-300, 400-410`

> `input: files array filenames including path of dG unfold files or ddG binding files used for calculation`

> `output: file including matrix of binding energies, or coupling energies`
'''

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'ENERGY', description = 'Generates a matrix of folding energies based on running energy minimization for mutational scanning.')
    parser.add_argument('--filenames', nargs = '+', type = str, dest = 'filenames', required = True, help = 'Filenames including path of dG unfold files or ddG binding files used for calculation', metavar = 'FILENS')
    parser.add_argument('--binding', type=str, required=True, dest = 'binding', help='String value to calculate binding ddG', metavar = 'BIND')
    parser.add_argument('--coupling', type=str, required=True, dest = 'coupling', help='String value to calculate coupling dddGs', metavar = 'COUPL')
    parser.add_argument('--llist_filename', type=str, default=None, help = 'Filename of location ranges of pdb locations to constrain energy calculations, comma delimited in file', metavar = 'LFILEN')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser


def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/energy/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')



''' Private functions '''
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



'''Private'''
def constrain_estimator_features(data, cutoffmin =0, cutoffmax=0, topmuts = 10, filter_name='chain', filter='H'):

    """ Constrain the estimator features 
    """

    indices_wt = (data['wild_type'].astype(bool) == True) & (data['coefficient'].astype(float) >=cutoffmax)
    indices_not_wt = (data['wild_type'].astype(bool) == False) & (data['coefficient'].astype(float) <= cutoffmin)
    indices_filter = (data[filter_name] == filter)
    indices = indices_wt + indices_not_wt
    mut= data.loc[indices]
    mut = mut.loc[indices_filter]
    mut = mut.loc[mut['coefficient'].astype(float).nlargest(n=topmuts).index]
    mut['location'] = mut['location'].astype(str)

    pdb, locations = read_pdb_locations(antibody=data['antibody_id'].values[0])

    merged = mut.merge(locations, how='left', left_on =['chain', 'location'], right_on = ['chain', 'fasta_location'])
    merged['mut'] = merged.wt_pdb + merged.aa_x

    return merged


'''Private'''
def constrain_energies(data, cutoffmin =0, cutoffmax=0, topmuts=10):
    

    data['mutation'] = data.mut.str[-1:]
    data = data.loc[~data.mutation.isin(['e','o'])]
    sorted = data.sort_values(by='ddg_bind').reset_index()
    constrained = sorted.loc[sorted.ddg_bind < cutoffmax ]
    constrained['loc'] = constrained['mut'].str[4:-1]
    constrained = constrained.iloc[0:topmuts+1,:]

    return constrained


'''Private'''    
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
    df['aa_three'] = [utils.aa_1to3_dic[aa] for aa in df.aa.values]
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


''' Private ''' 
def build_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4):

    mutations = get_optimal_antibody_structure_mutations(ab_name, p, p_less, upper_bound_ddg = -0.4)
    write_to_mutations_file(ab_name, mutations, 'individual_list.txt')
    run_build_model(ab_name, p + '.pdb', 'individual_list.txt')
    run_build_model(ab_name, p_less + '.pdb', 'individual_list.txt')

    run_multiple_repair_model(ab_name,p, mutations)
    run_multiple_repair_model(ab_name,p_less, mutations)
    return mutations


''' Private '''
def write_to_mutations_file(mutations, mutations_filename, out_dir):


    ab_path = '../../data/pdb/' + ab_name + '/'
    out_path = '../../data/pdb'  + ab_name + '/'

    mutstr = ''
    for mut in mutations:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutations_filename, 'w')   
    f.write(mutstr)
    f.close()

    f  = open(out_dir+ 'mutations.txt', 'w')
    f.write(mutstr)
    f.close()

    
'''API'''
def calculate_ddg_bind(antibody, pdb, pdb_less, scantype='virus', indir='./', outdir=''):
    ''' Calculate difference in Gibbs energy of binding for mutational scanning on virus bound to antibody
        antibody: name of antibody 
        pdb: PDB name, complexed structure 
        pdb_less:  PDB structure name, single uncomplexed structure
        scantype:  whether virus scanning or antibody scanning was performed 
        indir: input directory
        output: ouptut directory 
    '''

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
    ddg['mut'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[-4:]

    ddg['wt'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:4] + ddg.pdb_location
    ddg['chain'] = ddg.mut_chain.str[3]
    ddg['wildtype'] = ddg.mut.str[0] == ddg.mut.str[-1]

    ddg_name = 'ddg_bind'
    ddg[ddg_name] = ddg['dg'] - ddg['dg_less']
    ddg['antibody'] = antibody
    
    ''' Output the file to the appropriate directory ''' 
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ddg.to_csv(outdir + 'ddg_bind_' + antibody + scantype + '_scanning.csv', index=None)
    print (ddg.head())
    return ddg


''' API ''' 
def calculate_dddg_bind(antibody, pdb_opt, pdb_opt_less, pdb, pdb_less, scantype='virus', indir='./', outdir=''):
    ''' Calculate difference in Gibbs energy of binding for virus mutational scanning bound to antibody
        and designed antibody.  Differences in binding energy are calculated for places where designed
        mutations are <=4A from virus receptor.  
        antibody: name of antibody 
        pdb_opt: PDB name, complexed structure of optimized/designed antibody 
        pdb_opt_less:  PDB structure name of optimized/designed antibody, single uncomplexed structure 
        pdb: PDB name, complexed structure 
        pdb_less:  PDB structure name, single uncomplexed structure
        scantype:  whether virus scanning or antibody scanning was performed 
        indir: input directory
        output: ouptut directory 
    '''

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
    ddg['mut'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[-4:]

    ddg['wt'] = su.convert_pdb_to_fasta_style(ddg.mut_chain.str[0:3]) + ddg.mut_chain.str[3:4] + ddg.pdb_location
    ddg['chain'] = ddg.mut_chain.str[3]
    ddg['wildtype'] = ddg.mut.str[0] == ddg.mut.str[-1]

    ddg_name = 'ddg_bind'
    ddg[ddg_name] = ddg['dg'] - ddg['dg_less']
    
    ''' Output the file to the appropriate directory ''' 
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ddg.to_csv(outdir + 'ddg_bind_' + antibody + scantype + '_scanning.csv', index=None)
    print (ddg.head())
    return ddg

def check_args(namespace):
    # TODO
    pass
    

def run_energy(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)
    
    ###### IMPLEMENTATION #######
    
    # `input: filenames including path of dG unfold files or ddG binding files used for calculation, comma delimited` ===> namespace.filenames (this is an array of strings which are each filename)
    
    # `input: locations list of locations to constrain energy calculations, comma delimited, eg: 250-300, 400-410 ===> namespace.llist_filename, the string filename of list of range of locations

    #  `input: binding str to calculate binding ddG ` ===> namespace.binding

    #  `input: coupling str to calculate coupling dddGs ` ===> namespace.coupling
    # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/energy/) and figs in namespace.output_dir/output/figs/
