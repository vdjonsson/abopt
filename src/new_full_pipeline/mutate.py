import energy as e
import foldx as foldx
import pandas as pd
import structure
import os

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'MUTATE', description = 'Mutate an antibody given a list of mutations and generate a structure for each mutation.')
    parser.add_argument('--pdb', type = str, required = True, dest = 'pdb_directory', help = 'PDB number ', metavar = 'PDB')
    parser.add_argument('--pdb_directory', type = str, required = True, dest = 'pdb_directory', help = 'PDB directory', metavar = 'PDB_DIR')
    parser.add_argument('--o', type = str, required = True, dest = 'output_directory', help = 'Output directory', metavar = 'OUT')
    return parser

def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/mutate/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def mutate (pdb, mutations, pdb_dir, out_dir):

    foldx.create_individual_list(mutations,pdb,out_dir)
    foldx.run_build_model(pdb, 'individual_list.txt', pdb_dir, out_dir)
    foldx.rename_buildmodel_files(pdb[:-4], out_dir, './individual_list.txt')

def run_mutate(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    create_output_structure(namespace.output_directory)
    mutations = [] #TODO ??
    mutate(namespace.pdb, mutations, namespace.pdb_directory, namespace.output_directory+'output/mutate/')