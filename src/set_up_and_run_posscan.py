#!/usr/bin/env python
# Tea Freedman-Susskind
# Jonsson Computational Lab
# 28 July 2020


import os
import sys

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))

pdb_files = os.getcwd()

"""
Format to Call:
python3 set_up_and_run_posscan.py <pdb name> <chain name>

Example:
python3 set_up_and_run_posscan.py 6xcm H 
"""

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


def run_pos_scan(pdb_name, mut_list, repair = False):
    os.system("cd " + pdb_files)
    big_string = ""
    for item in mut_list:
        big_string += item + ","
    big_string = big_string[:-1]

    if repair: 
        repair = "foldx --command=RepairPDB --pdb="+pdb_name+".pdb"
        os.system(repair)

    command = "foldx --command=PositionScan --pdb="+pdb_name+"_Repair.pdb --positions="+big_string +" --out-pdb=false"
    os.system("echo " + big_string)
    os.system(command)



if __name__ == "__main__":
    args = sys.argv
    pdb_len = 1
    chain_len = 2
    repair = 3 
    if len(args) > pdb_len:
        pdb_name = args[pdb_len]
        file_name = pdb_name + '.pdb'
        os.system("echo " + pdb_name)
        l = pdb_files + file_name
        if len(args) > chain_len:
                chain_name = sys.argv[chain_len]
                run_pos_scan(pdb_name, pdb_to_list(file_name, spec_chain=chain_name))
        else:
                run_pos_scan(pdb_name, pdb_to_list(file_name))

