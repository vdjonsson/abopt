import os
import pandas as pd
#import wget
from biopandas.pdb import PandasPdb


aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL','NA']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V' ,'X']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))

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



def delete_specified_chains(atom_panda, labeled_chains, keyword_list):# atom_panda, labeled_chains, keyword_list):
    """
    Remove given chains from a given pdb biopanda
    :param atom_panda: a biopanda containing pdb data
    :param labeled_chains: dictionary matching chain names to pdb tags
    :param keyword_list: words indicating a chain to delete
    :return: biopanda with desired chains removed
    """
    keyword_list = [value.lower() for value in keyword_list]
    chains_to_delete = [value for key, value in labeled_chains.items() if any(item in key.lower() for item in keyword_list)]
    flat_chains = [item for sublist in chains_to_delete for item in sublist]
    reduced = atom_panda[~atom_panda.chain_id.isin(flat_chains)]
    return reduced


def remove_virus_with_biopandas(pdb_file_name, labeled_chains):
    """
    Removes SARS-CoV-2 spike glycoprotein chains from a pdb and prints to a new pdb file
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :return: biopanda containing pdb without the virus
    """
    struc_path = os.path.abspath( pdb_file_name + ".pdb")
    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    the_pdb = ppdb.df['ATOM']
    the_pdb = delete_specified_chains(the_pdb, labeled_chains, ['spike'])
    ppdb.df['ATOM'] = the_pdb
    file_path = pdb_file_name + '_no_virus.pdb'
    ppdb.to_pdb(path=file_path, records=['ATOM'], gz=False, append_newline=True)
    return the_pdb


def remove_chains(pdb_file_name, labeled_chains, to_del):
    """
    Removes input chains from a pdb and prints to a new pdb file
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :param to_del: list of keywords indicating chains to delete
    :return: biopanda containing pdb without the chains
    """
    struc_path = os.path.abspath(data_path + pdb_file_name + ".pdb")
    ppdb = PandasPdb()
    ppdb.read_pdb(struc_path)
    the_pdb = ppdb.df['ATOM']
    res = delete_specified_chains(the_pdb, labeled_chains, to_del)
    ppdb.df['ATOM'] = res
    file_path = data_path + pdb_file_name + '_reduced.pdb'
    ppdb.to_pdb(path=file_path, records=['ATOM'], gz=False, append_newline=True)
    return ppdb.df['ATOM']


def find_epitopes(pdb_file_name, labeled_chains):
    """
    Find epitopes (ab to rbd <= 4A) on the antibody, output to csv
    :param pdb_file_name: name of pdb
    :param labeled_chains: dictionary matching chain names to pdb tags
    :return: panda containing antibody epitope locations
    """
    struc_path = os.path.abspath(data_path + pdb_file_name + ".pdb")
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
        spike_within_4A = ppdb.df['ATOM'][distances < 4.0]
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
    out_file = data_path + pdb_file_name + '_epitopes.csv'
    epitopes.to_csv(out_file, index=False)
    return epitopes


def convert_pdb_to_fasta_style(pdb_mutations):
    mut_one = []
    for mut in pdb_mutations:
        if mut !='H2S':
            mut_one.append(aa_3to1_dic[mut])
    return pd.Series(mut_one)

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


def test(inpt):
    mut_panda = pd.read_csv(inpt)
    mut_dict = {}
    for index, row in mut_panda.iterrows():
        if not row['wildtype']:
            pdb_name = row['pdb']
            mut = row['mut']
            if pdb_name in mut_dict:
                mut_dict[pdb_name].append(mut)
            else:
                mut_dict[pdb_name] = [mut]
    for pdb in mut_dict:
        temp = mut_panda.loc[mut_panda['pdb'] == pdb]
        temp_ab = temp.iloc[0, :]
        ab = temp_ab['antibody']
        indiv_list = data_path + 'individual_list_' + pdb + '.txt'
        write_to_mutations_file(mut_dict[pdb], indiv_list)
        pdb_url = 'https://files.rcsb.org/download/' + pdb + '.pdb'
        pdb_file = wget.download(pdb_url)
        run_repair_model(ab, pdb.lower())
        rep_pdb = pdb.lower() + '_Repair'
        run_build_model(ab, rep_pdb, indiv_list)
        rename_bm_out(rep_pdb, indiv_list)
        for mut in mut_dict[pdb]:
            mut_pdb = rep_pdb + '_' + mut
            print(mut_pdb)
            print(label_chains(pdb.lower()))
            rm_virus_with_biopandas(mut_pdb, label_chains(pdb.lower()))
            no_vir_pdb = mut_pdb + '_no_virus'
            run_repair_model(ab, no_vir_pdb)

# inp = '/home/teafs/Downloads/tmp.csv'
# test(inp)
