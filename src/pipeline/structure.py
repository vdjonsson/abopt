import os
import pandas as pd
import wget
from biopandas.pdb import PandasPdb
import utils as utils


def setup_parser():
    parser = argparse.ArgumentParser(prog = 'STRUCTURE', description = 'Functions to analyze PDB structures') 
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of fitness matrix dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of fitness matrix dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--p', type = float, required = True, dest = 'p_unmanaged', help = 'Maximum proportion of viruses not covered by the antibody cocktail', metavar = 'P_UNMANAGED')
    parser.add_argument('--l', type = float, required = True, dest = 'lmbda', help = 'Lambda value to use for penalization of number of antibodies chosen', metavar = 'LMBD')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser


def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+ 'output/structure/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')





''' Private ''' 
def delete_specified_chains(atom_panda, labeled_chains, keyword_list):# atom_panda, labeled_chains, keyword_list):
    """
    Remove given chains from a given pdb biopanda
    :param atom_panda: a biopanda containing pdb data
    :param labeled_chains: dictionary matching chain names to pdb tags
    :param keyword_list: words indicating a chain to delete
    :return: biopanda with desired chains removed
    """
    keyword_list = [value.lower() for value in keyword_list]
    print('keyword list')
    print(keyword_list)
    chains_to_delete = [value for key, value in labeled_chains.items() if any(item in key.lower() for item in keyword_list)]
    print('chains to delete')
    print(chains_to_delete)
    flat_chains = [item for sublist in chains_to_delete for item in sublist]
    print(flat_chains)
    reduced = atom_panda[~atom_panda.chain_id.isin(flat_chains)]
    return reduced


''' API ''' 
def label_chains(pdb_name):
    """
    Match chain names to alphabetical tags in the PDB file
    :param pdb_name: the name of the file on the RCSB PDB library
    :return: dictionary matching chain names to pdb tags
    """
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



''' API ''' 
def remove_virus_with_biopandas(pdb_file_name, labeled_chains):
    """
    Removes SARS-CoV-2 spike glycoprotein chains from a pdb and prints to a new pdb file
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :return: biopanda containing pdb without the virus
    """
    struc_path = os.path.abspath( pdb_file_name)
    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    the_pdb = ppdb.df['ATOM']
    the_pdb = delete_specified_chains(the_pdb, labeled_chains, ['spike'])
    ppdb.df['ATOM'] = the_pdb
    file_path = pdb_file_name + '_no_virus.pdb'
    ppdb.to_pdb(path=file_path, records=['ATOM'], gz=False, append_newline=True)
    return the_pdb


''' API''' 
def remove_chains(pdb_dir,pdb, labeled_chains, chain_type, out_dir):
    """
    Removes input chains from a pdb and prints to a new pdb file
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :param to_del: list of keywords indicating chains to delete
    :return: biopanda containing pdb without the chains
    """
    less_molecule = 'less_virus'
    chain_keyword = ['Spike']
    if chain_type == 'antibody': 
        less_molecule = 'less_ab'
        chain_keyword = ['Fab', 'chain']
    struc_path = os.path.abspath(pdb_dir + pdb )

    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    the_pdb = ppdb.df['ATOM']
    res = delete_specified_chains(the_pdb, labeled_chains, chain_keyword)

    ppdb.df['ATOM'] = res
    file_path = out_dir + pdb[:-4] + '_' + less_molecule + '.pdb'
    print('Removed: ' + file_path)
    ppdb.to_pdb(path=file_path, records=['ATOM'], gz=False, append_newline=True)
    return ppdb.df['ATOM']


def find_epitopes(pdb_dir, pdb_file_name, labeled_chains, distance):
    """
    Find epitopes (ab to rbd <= distance) on the antibody, output to csv
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :return: panda containing antibody epitope locations
    """
    struc_path = os.path.abspath(pdb_dir + pdb_file_name)
    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    ab_lis = [key for key, value in labeled_chains.items() if 'spike' not in key.lower()]
    spike_lis = [key for key, value in labeled_chains.items() if 'spike' in key.lower()]
    spike_chains = [value for key, value in labeled_chains.items() if 'spike' in key.lower()]
    ab = delete_specified_chains(ppdb.df['ATOM'], labeled_chains, spike_lis)
    ppdb.df['ATOM'] = delete_specified_chains(ppdb.df['ATOM'], labeled_chains, ab_lis)
    cols = ['chain_antibody', 'number_antibody', 'residue_antibody', 'closest distance', 'chain_virus', 'number_virus', 'residue_virus']
    epitopes = pd.DataFrame(columns=cols)
    i = 0
    prev_res = ''
    for index, row in ab.iterrows():
        ref_pt = (row['x_coord'], row['y_coord'], row['z_coord'])
        distances = ppdb.distance(xyz=ref_pt, records=('ATOM',))
        spike_within_4A = ppdb.df['ATOM'][distances < distance]
        res = row['residue_number']
        if not spike_within_4A.empty and res != prev_res:
            min_ind = distances.idxmin()
            minm = distances[min_ind]
            closest_spike = spike_within_4A.loc[min_ind, :]
            new_row = {'chain_antibody': row['chain_id'], 'number_antibody': res,
                       'residue_antibody': row['residue_name'], 'closest distance': minm,
                       'chain_virus': closest_spike['chain_id'], 'number_virus': closest_spike['residue_number'],
                       'residue_virus': closest_spike['residue_name']}
            epitopes.loc[i] = new_row
            i += 1
            prev_res = res
    out_file = '../../output/epitope/' + pdb_file_name[:-4] + '_epitopes.csv'
    epitopes.to_csv(out_file, index=False)
    return epitopes


''' API '''
def convert_pdb_to_fasta_style(pdb_mutations):
    mut_one = []
    for mut in pdb_mutations:
        mut_one.append(utils.aa_3to1_dic[mut])
    return pd.Series(mut_one)

''' API '''
def convert_pdb_to_list(ab_name, pdb_name, spec_chain='C' , min_res=400, max_res=520):
    ''' Convert PDB location to a list of mutations 
    '''

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
                    for amac in utils.aa_3to1_dic:
                        if amac != broken[aa]:
                            try:
                                utils.aa_3to1_dic[broken[aa]]
                            except KeyError:
                                print(broken)
                                for amac in utils.aa_3to1_dic:
                                    if amac in broken[aa]:
                                        print(amac)
                                        broken[aa] = amac
                                        print(broken)
                        to_list = utils.aa_3to1_dic[broken[aa]] + broken[chain] + broken[residue_no] + 'a'

                        if (int(broken[residue_no]) >= min_res and int(broken[residue_no])<= max_res):
                            lis.append(to_list)

                        prevnum = broken[residue_no]
    return lis

