#!/usr/bin/env python
# Tea Freedman-Susskind
# Jonsson Computational Lab
# 16 July 2020


import os

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))

pdb_files = '/path/to/pdb/directory'


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
                    lis.append(to_list)
                    prevnum = broken[residue_no]
    return lis


def run_pos_scan(pdb_name, mut_list):
    os.system("echo what up")
    os.system("cd /path/to/pdb/directory")
    big_string = ""
    for item in mut_list:
        big_string += item + ","
    big_string = big_string[:-1]
    command = "foldx --command=PositionScan --pdb="+pdb_name+".pdb --positions="+big_string +" --out-pdb=false"
    os.system("echo " + big_string)
    os.system(command)


l = pdb_files + 'pdbname.pdb'
run_pos_scan('pdbname', pdb_to_list(l))
