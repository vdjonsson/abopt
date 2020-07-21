#!/usr/bin/env python
# Tea Freedman-Susskind
# Jonsson Computational Lab
# 21 July 2020

import seaborn as sns
import matplotlib.pyplot as plt
import csv
import pandas as pds

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))

pdb_files = '/path/to/pdb/files'


def pdb_to_output(pdb, outname, spec_chain=''):
    indivlist = pdb_files + 'individual_list' + outname + '.txt'
    prevnum = 0
    type_id = 0
    aa = 3
    chain = 4
    residue_no = 5
    with open(pdb, 'r') as f:
        with open(indivlist, 'w') as g:
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
                                to_file = aa_3to1_dic[broken[aa]] + broken[chain] + broken[residue_no] + aa_3to1_dic[amac] + ";\n"
                                g.write(to_file)
                    prevnum = broken[residue_no]


def ddg_from_output(pdb):
    filename = pdb_files + 'PS_' + pdb + '_scanning_output.txt'
    out_dict = {}
    mutname = 0
    wt = 3
    with open(filename, 'r') as f:
        for line in f:
            broken = line.split("	")
            wt_aa = broken[mutname][:wt]
            if wt_aa.startswith('H') and wt_aa.endswith('S'):
                wt_aa = 'HIS'
            key = aa_3to1_dic[wt_aa] + broken[mutname][wt:]
            out_dict[key] = float(broken[1].strip())
    return out_dict


def make_ddg_csv(pdb, ddg_dict):
    new_csv = pdb_files + pdb + ".csv"
    with open(new_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['wt/chain/location/mutation','ddG'])
        for muta in ddg_dict:
            writer.writerow([muta, ddg_dict[muta]])


def chain_sep_list(ddgs):
    list_of_lists = []
    current_chain = 0
    ind = -1
    chain_ind = 1
    for mut in ddgs:
        chain = mut[chain_ind]
        if chain == current_chain:
            list_of_lists[ind].append([mut, ddgs[mut]])
        else:
            ind += 1
            list_of_lists.append([[mut, ddgs[mut]]])
            current_chain = chain
    return list_of_lists


def chain_sep_dict(ddgs):
    list_of_dicts = []
    current_chain = 0
    ind = -1
    chain_ind = 1
    for mut in ddgs:
        chain = mut[chain_ind]
        if chain == current_chain:
            list_of_dicts[ind][mut] = ddgs[mut]
        else:
            ind += 1
            list_of_dicts.append({mut: ddgs[mut]})
            current_chain = chain
    return list_of_dicts


def pdb_to_list(pdb, spec_chain=''):
    prevnum = 0
    type_id = 0
    aa = 3
    chain = 4
    residue_no = 5
    lis = []
    with open(pdb, 'r') as f:
        for line in f:
            broken = line.split(" ")
            broken = list(filter(None, broken))
            if broken[type_id] == "ATOM" and broken[residue_no] != prevnum and broken[residue_no].isnumeric():
                if not spec_chain or (spec_chain and spec_chain == broken[chain]):
                    to_list = broken[aa] + broken[chain] + broken[residue_no] + 'a'
                    lis.append(to_list)
                prevnum = broken[residue_no]
    return lis

pdb_name = '6xdg_REGN10933_RBD_only'
ddg_dict = ddg_from_output(pdb_name)
make_ddg_csv(pdb_name, ddg_dict)
fig_dir = "/path/to/where/you/want/figures"
# """
# MAKE FREQUENCY PLOT OF FULL STRUCTURE DDGs OF MUTATION
ddg_list = list(ddg_dict.values())
ddg_freq = sns.distplot(ddg_list)
ddg_freq.set_title(pdb_name + " ddG mutation frequency")
ddg_freq.set(xlabel='ddG, kcal/mol', ylabel='freq')
plt.savefig(pdb_name + "_freq_all.png")
plt.close()
# """
chains_ddg = chain_sep_list(ddg_dict)
chains_ddg_dict = chain_sep_dict(ddg_dict)
for chain in chains_ddg:
    chain_pd = pd.DataFrame(chain, columns=['mutation', 'ddG(kcal/mol)'])
    chain_ind = 1
    chain_name = chain_pd['mutation'][0][chain_ind]
    # """
    # MAKE FREQUENCY PLOT OF EACH CHAIN DDGs OF MUTATION
    chain_freq = sns.distplot(chain_pd['ddG(kcal/mol)'])
    chain_freq.set_title(pdb_name + " " + chain_name + " chain \u0394\u0394G mutation frequency")
    chain_freq.set(xlabel='ddG, kcal/mol', ylabel='freq')
    plt.savefig(fig_dir + pdb_name + "_" + chain_name + "_ddg_freq.png")
    plt.close()
    # """

for chain in chains_ddg_dict:
    by_location = {}
    chain_ind = 1
    keycode = '3000'
    for mut in chain:
        if keycode in mut:
            by_location[keycode].append(chain[mut])
        else:
            keycode = mut[:-1]
            by_location[keycode] = [chain[mut]]
    chain_name = keycode[chain_ind]
    # """
    # MAKE BOXPLOT OF DDGS AT EACH LOCATION
    by_loc_pd = pd.DataFrame.from_dict(by_location, orient='index')
    by_loc_pd = by_loc_pd.transpose()
    plt.figure(figsize=(20, 6.5))
    by_loc_box = sns.boxplot(data=by_loc_pd)
    by_loc_box.set_title(pdb_name + " " + chain_name + " chain \u0394\u0394G Mutation Distribution by Location")
    by_loc_box.set(xlabel='residue', ylabel='ddG, kcal/mol')
    for tick in by_loc_box.xaxis.get_major_ticks():
        tick.label.set_fontsize(5)
    plt.tight_layout()
    plt.xticks(rotation=90)
    plt.savefig(fig_dir + pdb_name + "_" + chain_name + "_ddg_box.png")
    plt.close()
    # """

