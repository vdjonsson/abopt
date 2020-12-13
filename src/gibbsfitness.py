#!/usr/bin/env python
# Tea Freedman-Susskind
# Jonsson Computational Lab
# 30 June 2020
# DEPRECATED 

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_dict = dict(zip(aa_single, aa_three))
aa_tcid = dict(zip(aa_three, aa_single))

pdb_files = '/path/to/pdb/files'


def pdb_to_output(pdb, outname):
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
                    for amac in aa_tcid:
                        if amac != broken[aa]:
                            try:
                                aa_tcid[broken[aa]]
                            except KeyError:
                                print(broken)
                                for amac in aa_tcid:
                                    if amac in broken[aa]:
                                        print(amac)
                                        broken[aa] = amac
                                        print(broken)
                            to_file = aa_tcid[broken[aa]] + broken[chain] + broken[residue_no] + aa_tcid[amac] + ";\n"
                            g.write(to_file)
                    prevnum = broken[residue_no]


def return_matched_ddg_from_fxout(fxout, mutant_file):
    dict = {}
    with open(fxout, 'r') as f:
        with open(mutant_file, 'r') as g:
            ddg = 2
            for line in f:
                broken = line.split("\t")
                if len(broken) > 3 and 'total energy' not in broken[ddg]:
                    item_ddg = broken[ddg]
                    mut = g.readline().strip("\n").strip(";")
                    dict[mut] = item_ddg
    return dict



pdb = '7bz5'
filename = pdb_files + pdb + '.pdb'
fil = pdb_files + "Average_" + pdb + ".fxout"
A = pdb_files + '7bz5_Chain_A.pdb'
pdb_to_output(A, 'chainA')
H = pdb_files + '7bz5_Chain_H.pdb'
pdb_to_output(H, 'chainH')
L = pdb_files + '7bz5_Chain_L.pdb'
pdb_to_output(L, 'chainL')
