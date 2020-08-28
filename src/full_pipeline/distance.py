import pandas as pd
import numpy as np
import argparse
from Levenshtein import distance as levenshtein_distance

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'DISTANCE', description = 'Calculate pairwise Levenshtein (edit) distance for provided sequences.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser
    
def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/distance/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def check_dataframe(df):
    if 'heavy_chain' not in df.columns:
        raise ValueError('Dataframe object must have column with heavy chain sequences named "heavy_chain".')
    if 'light_chain' not in df.columns:
        raise ValueError('Dataframe object must have column with light chain sequences named "light_chain".')
    if 'antibody_id' not in df.columns:
        raise ValueError('Dataframe object must have column identifying each sequence named "antibody_id".')

def calculate_distances(df, id_col = 'antibody_id', heavy_col = 'heavy_chain', light_col = 'light_chain'):
    l_distances = pd.DataFrame(data = np.array([np.nan]*(len(df)**2)).reshape((len(df), len(df))), columns = df[id_col].values, index = df[id_col].values)

    col_combns = combinations(df[id_col].values, 2)
    for cols in list(col_combns):
        seq1_heavy, seq1_light = df.loc[df[id_col] == cols[0]][[heavy_col, light_col]].values.flatten()
        seq2_heavy, seq2_light = df.loc[df[id_col] == cols[1]][[heavy_col, light_col]].values.flatten()
        l_dist = levenshtein_distance(seq1_heavy, seq2_heavy) + levenshtein_distance(seq1_light, seq2_light)
        l_distances[cols[0]][cols[1]] = l_dist
        l_distances[cols[1]][cols[0]] = l_dist

    l_distances[np.isnan(l_distances)] = 0
    return l_distances

def run_distance(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    create_output_structure(namespace.output_dir)
    
    df = pd.read_csv(namespace.filepath+namespace.filename+'.csv', sep=',', header=0)
    check_dataframe(df)
    l_distances = calculate_distances(df)
    
    l_distances.to_csv(namespace.output_dir+'output/distance/'+namespace.filename+'.csv', sep=',', header=True, index=True)
