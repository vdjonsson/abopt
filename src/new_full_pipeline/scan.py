import energy as e
import foldx as foldx
import pandas as pd
import structure
import os 

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'SCAN', description = 'Mutational scanning on a subset of locations or a chain.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of dataframe', metavar = 'FILEP')
    # TODO
    return parser

def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/scan/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def scan (scantype, scanvalues, scanmolecule, antibody, pdblist, pdbdir, outdir):

    """ Mutational scanning
    """

    if scantype == 'mutation':
        pdbname, pdb_loc = e.read_pdb_locations(file_location=file_ab_locations)
    elif scantype =='location':
        posscan_str =  ",".join(scanvalues)
    elif scantype == 'chain': # position scan entire chain
        posscan_str =  ","

    for pdb in pdblist:
        foldx.run_position_scan (pdb, scanmolecule,  posscan_str, pdbdir, outdir)

def run_scan(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    create_output_structure(namespace.output_directory)
    # TO DO
