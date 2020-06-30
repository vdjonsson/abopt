#!/usr/bin/env python
# Tea Freedman-Susskind
# Jonsson Computational Lab
# 30 June 2020

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_dict = dict(zip(aa_single, aa_three))
aa_tcid = dict(zip(aa_three, aa_single))

pdb_files = '/directory/pdb/files/are/in'

def pdb_to_output(pdb):
    indivlist = pdb_files + 'individual_list.txt'
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
                if broken[type_id] == "ATOM" and broken[residue_no] != prevnum and len(broken) > 12:
                    for amac in aa_tcid:
                        if amac != broken[aa]:
                            try:
                                aa_tcid[broken[aa]]
                            except KeyError:
                                for amac in aa_tcid:
                                    if amac in broken[aa]:
                                        broken[aa] = amac
                            to_file = aa_tcid[broken[aa]] + broken[chain] + broken[residue_no] + aa_tcid[amac] + ";\n"
                            g.write(to_file)
                prevnum = broken[residue_no]


pdb = '7bz5'
filename = pdb_files + pdb + '.pdb'
pdb_to_output(filename)
