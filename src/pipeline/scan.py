import foldx as fx
import argparse

'''

> `input: foldx use foldx for mutational scanning `
> `input: filenames of molecular structures to scan, comma delimited`
> `input: scantype virus or antibody`
> `input: chain scan the entire chain`
> `input: locations range of pdb locations to scan, comma delimited, eg: 400-403,420-423`
> `output: file with dG of unfolding from molecular structure`
'''

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'SCAN', description = 'Performs mutational scanning on a subset of locations or a chain.')
    parser.add_argument('--backend', type = str, required = False, dest = 'backend', default = 'foldx', help = 'Software program to use as backend for mutational scanning (i.e. FoldX, Rosetta)', metavar = 'PKG')
    parser.add_argument('--filenames', nargs = '+', type = str, dest = 'scantype', required = True, help = 'Filenames of molecular structures to scan', metavar = 'FILENS')
    parser.add_argument('--scantype', choices = ['antibody', 'virus'], type=str, required = True, dest = 'scantype', help = 'Type of structure to scan ("antibody" or "virus")', metavar = 'TYPE')
    parser.add_argument('--chain', action = 'store_const', const = True, default = False, dest = 'chain', help = 'Flag to indicate to scan the entire chain')
    parser.add_argument('--llist_filename', type=str, default=None, help = 'Filename of location ranges of pdb locations to scan, comma delimited in file', metavar = 'LFILEN')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
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
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')


def check_args(namespace):
    # TODO
    if namespace.backend.lower() != 'foldx':
        raise NotImplementedError
    pass

def run_scan(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)
    
    ###### IMPLEMENTATION #######
    
    # `input: foldx use foldx for mutational scanning ` ===> namespace.backend
    # `input: filenames of molecular structures to scan, comma delimited` ===> namespace.filenames (this is an array of strings which are each filename)
    # `input: scantype virus or antibody` ===> namespace.scantype
    # `input: chain scan the entire chain` ===> namespace.chain (boolean value)
    # `input: locations range of pdb locations to scan, comma delimited, eg: 400-403,420-423` ===> namespace.llist_filename, the string filename of list of range of locations to scan
    # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/scan/) and figs in namespace.output_dir/output/figs/

    if namespace.backend == 'foldx':

        for fname in filenames: 
            pdb_file = fname[fname.rfind('/'):]
            pdb_dir  = fname[:fname.rfind('/')]
            pos_scan_str = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)

            # assume you have a location file for this antibody 
            # pos_scan_str = fx.construct_positionscan_string (pdb_name, chain = namespace.chain, location_file)
            pos_scan_str = ''

            fx.run_position_scan (pdb_file=pdb_file, scan_molecule = namespace.scantype, pos_scan=pos_scan_str, pdb_dir = pdb_dir, out_dir = namespace.output_dir)


