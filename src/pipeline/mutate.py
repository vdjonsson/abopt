import foldx as fx
import argparse

''' 

abopt mutate mutates an a structure given a list of mutations, and generates a structure for each mutation/mutation list.

input: filename str of molecular structure to mutate

input: chain str of the chain to mutate used when mutating locations

input: mutations array list of mutations or locations comma delimited, eg: TH28I, YH58F

input: location array list of locations to mutate, all amino acids comma delimited, eg: 28-58, 75

input: repair bool True if structure(s) to mutate requires repair after mutating

output: PDB files of mutated and repaired molecular structures


'''

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'MUTATE', description = 'Mutate a molecular structure')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filepath/filename of molecular structures to mutate', metavar = 'FILEN')
    parser.add_argument('--chain', type = str, required = True, dest = 'chain', help = 'Chain on structure to mutate used when mutating locations')
    parser.add_argument('--mlist_filename', type = str, default = None, dest = 'mlist_filename', help = 'Filepath/filename for list of mutations', metavar ='MFILEN')
    parser.add_argument('--llist_filename', type = str, default = None, dest = 'llist_filename', help = 'Filepath/filename for list of locations to mutate', metavar ='LFILEN')
    parser.add_argument('--repair', action='store_const', const = True, default = False, dest = 'repair', help = 'Flag to use if structure(s) to mutate require(s) repair after mutation')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
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
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def run_mutate(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    create_output_structure(namespace.output_dir)
    
    ## IMPLEMENTATION HERE: where to find variables:
    
    # input: filename str of molecular structure to mutate ==> namespace.filename
    # input: chain str of the chain to mutate used when mutating locations ==> namespace.chain
    # input: mutations array list FILE of mutations or locations comma delimited, eg: TH28I, YH58F ==> namespace.mlist_filename
    # input: location array list FILE of locations to mutate, all amino acids comma delimited, eg: 28-58, 75 ==> namespace.llist_filename
    # input: repair bool True if structure(s) to mutate requires repair after mutating ==> namespace.repair (this is a boolean value)

    # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/mutate/) and figs in namespace.output_dir/output/figs/

    # parse out the pdb and the path 
    pdb_name = namespace.filename[namespace.filename.rfind('/'):]
    pdb_dir = namespace.filename[:namespace.filename.rfind('/')]

    if namespace.mlist_filename !=None:
        fx.run_build_model(pdb_name=pdb_name, mutations_file=namespace.mlist_filename, pdb_dir=pdb_dir, out_dir = namespace.output_dir)
    elif namespace.llist_filename != None:
        raise NotImplementedError
