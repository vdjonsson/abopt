import energy as e
import foldx as foldx
import pandas as pd
import structure
import os
import argparse

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'CONSTRAIN', description = 'Constrain estimator features.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--name', type = str, required = True, dest = 'name', help = 'Name of antibody in dataframe to map to PDB locations', metavar = 'NAME')
    parser.add_argument('--name', type = str, required = True, dest = 'name', help = 'Name of antibody in dataframe to map to PDB locations', metavar = 'NAME')
    parser.add_argument('--type', choices=['estimator', 'energy'], dest = 'type', help = 'Type of constrain ("energy" or "estimator")', metavar = 'TYPE')
    parser.add_argument('--cutoff', type = int, action='append', required = True, dest = 'cutoff', help = 'Minimum and maximum value for cutoff of coefficients', metavar = '[MIN MAX]')
    parser.add_argument('--top', type = int, required = True, dest = 'top', help = 'Top number of locations to consider', metavar = 'TOP')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser
    
def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/constrain/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def constrain(constraintype, constrainfile, antibody, cutoff, top, out_dir):

   """ Constrain data in constrainfile using cutoff
   """

   if constraintype == 'estimator':
       estimator = e.read_estimator(constrainfile, antibody)
       constrained = e.constrain_estimator_features(estimator, cutoffmin =cutoff[0], cutoffmax=cutoff[1], topmuts = top, filter_name='chain', filter='H')
       constrained.to_csv(out_dir + antibody + '_estimator.csv', index=None)
   elif constraintype == 'energy':
       energies = pd.read_table(constrainfile, sep=',')
       constrained = e.constrain_energies(energies, cutoffmin= cutoff[0], cutoffmax= cutoff[1], topmuts=top)
       constrained.to_csv(out_dir + antibody + '_energies.csv', index=None)

def run_constrain(args):
    parser = setup_parser()
    namespace = parser.parse_args()
    create_output_structure(namespace.output_dir)
    
    constrain(namespace.type, namespace.filepath+namespace.filename, namespace.name, namespace.cutoff, namespace.top, namespace.output_dir)
