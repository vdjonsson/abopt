import pandas as pd
import numpy as np
import argparse
from estimator_utils.io_utils import read_file

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'MAP', description = 'Map locations indexed from 0 to provided PDB locations')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--pdb_filepath', type = str, required = True, dest = 'pdb_filepath', help = 'Filepath of dataframe mapping FASTA locations to PDB locations', metavar = 'FILEP')
    parser.add_argument('--pdb_filename', type = str, required = True, dest = 'pdb_filename', help = 'Filename of dataframe mapping FASTA locations to PDB locations (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--name', type = str, required = True, dest = 'name', help = 'Name of antibody in dataframe to map to PDB locations', metavar = 'NAME')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser

def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/map/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def check_dataframe(df):
    if 'antibody_id' not in df.columns:
        raise ValueError('Dataframe object must have column identifying each antibody named "antibody_id."')
    if 'location' not in df.columns:
        raise ValueError('Dataframe object must have column identifying location of each amino acid named "location."')
    if 'chain' not in df.columns:
        raise ValueError('Dataframe object must have column identifying chain of each amino acid named "chain."')
    if 'coefficient' not in df.columns:
        raise ValueError('Dataframe object must have column identifying coefficient associated with each location named "coefficient."')

def check_pdb_dataframe(pdb_df):
    if 'pdb_location' not in df.columns:
        raise ValueError('Dataframe object must have column identifying PDB location of each amino acid named "pdb_location."')
    if 'fasta_location' not in df.columns:
        raise ValueError('Dataframe object must have column identifying FASTA location of each amino acid named "fasta_location."')
    if 'chain' not in df.columns:
        raise ValueError('Dataframe object must have column identifying chain of each amino acid named "chain."')

def run_map(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    create_output_structure(namespace.output_dir)
    
    df = read_file(namespace.filepath, namespace.filename)
    pdb_df = read_file(namespace.pdb_filepath, namespace.pdb_filename)
    check_dataframe(df)
    check_pdb_dataframe(pdb_df)
    
    antibody_df = df.loc[df.antibody_id == namespace.name]
    chain_locations = []
    for chain in np.unique(antibody_df.chain.values):
        location_dict = dict(zip(pdb_df.loc[pdb_df.chain==chain].fasta_location.values, pdb_df.loc[pdb_df.chain==chain].pdb_location.values))
        chain_locations.append(antibody_df.loc[antibody_df.chain == chain].replace({'location': location_dict}))

    pdb_antibody_df = pd.concat(chain_locations)
    
    pdb_antibody_df.to_csv(namespace.filename+'_'+namespace.name+'_pdb_locations.csv', sep=',', header=True, index=False)
