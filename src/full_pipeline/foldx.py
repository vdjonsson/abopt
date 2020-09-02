import os 
import pandas as pd 


foldx_path = '/Applications/foldx5MacC11/'

def construct_positionscan_string (pdb_name, ab_pos, pdb_loc):
    """ Construct a position scan string for FoldX input 
    """ 
    check = (ab_pos.antibody_id.unique() == pdb_loc.antibody_id.unique())
    if check == False:
        print('Not matching antibody')
        exit()

    merged = ab_pos.join(pdb_loc, on='fasta_location', how = 'left', lsuffix='_abpos')
    merged = merged.dropna(axis=0, subset=['pdb_location'])
    posscan = merged.pdb_location.unique()
    pd.Series(posscan).to_csv(data_path + pdb_name + '_position_scan.csv', index=False)
    posscan_str =  ",".join(posscan)
    return posscan_str


def run_position_scan (pdb_file, scan_molecule, pos_scan, pdb_dir, out_dir):

    """ Run position scanning using FoldX
    """ 
    tag = pdb_file[:-4] + '_' + scan_molecule

    print(tag)
    print('Running position scan')
    command = foldx_path + "foldx --command=PositionScan --pdb-dir=" + pdb_dir + " --pdb=" + pdb_file
    command = command + " --positions="+ pos_scan +" --out-pdb=false  --output-dir=" + out_dir 
    command = command + " --output-file=" + tag
    print(command)
    os.system(command)

def run_build_model(pdb_name, mutations_file, pdb_dir, out_dir):
    """ Mutate a structure given a mutations file by running BuildModel using FoldX 
    """
    command = foldx_path + "foldx --command=BuildModel --pdb-dir=" + pdb_dir + " --pdb=" + pdb_name
    command = command + " --mutant-file="+ mutations_file + " --output-dir=" + out_dir
    print(command)
    os.system(command)

def rename_pdb_files(pdb_name, mutations, pdb_dir, out_dir):
    """ Rename PDB files after FoldX 
    """
    print(pdb_name)
    i = 1 
    for mut in mutations:
        pdb =  out_dir 
        command  = 'mv ' + pdb + '_' + str(i) + '.pdb ' + pdb + '_' + mut + '.pdb'
        print(command)
        os.system(command)
        i = i+1



def create_individual_list(mutations, mutpath):
    mutstr = ''
    for mut in mutations:
        mutstr = mutstr + mut+ ';\n'            
    f  = open(mutpath + 'individual_list.txt', 'w')   
    f.write(mutstr)
    f.close()

    return mutstr


def rename_buildmodel_files(pdb_name, pdb_dir, indiv_list_path):
    """
    Renames files output from BuildModel to label with mutation
    :param pdb_name: original pdb that was mutated
    :param indiv_list_path: the mutation list used in BuildModel
    :return: None
    """
    with open(indiv_list_path, 'r') as f:
        full = f.read()
    broken = full.split(';\n')
    for ind in range(len(broken)):
        if broken[ind]:
            file_to_find = pdb_dir + pdb_name + '_' + str(ind+1) + '.pdb'
            new_name = pdb_dir + pdb_name + '_' + broken[ind] + '.pdb'
            os.rename(file_to_find, new_name)


def run_repair_model(pdb_name, pdb_dir, out_dir):

    """ Repair a structural model by running RepairPDB using FoldX
    """

    command = foldx_path + "foldx --command=RepairPDB --pdb-dir=" + pdb_dir + " --pdb=" + pdb_name
    command = command +  " --output-dir=" + out_dir
    print(command)
    os.system(command)

def run_multiple_repair_model(pdb_name, mutations, pdb_dir, out_dir):

    """ Repair multiple mutated structural models given a list of mutations by running RepairPDB using FoldX
    """
    i = 1
    for mut in mutations:
        print(mut)
        run_repair_model(pdb_name + '_' + mut +'.pdb',  pdb_dir, out_dir)
        i = i+1
