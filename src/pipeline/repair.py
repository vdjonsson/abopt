import foldx as fx
import argparse


'''
`abopt repair` repairs an antibody for downstream analysis. 

> `input: array of filenames of molecular structures to repair`
> `input: foldx use foldx `
> `input: path pathname of tool used for repair `
> `output: PDB file of molecular structure mutated`
'''


def setup_parser():
    parser = argparse.ArgumentParser(prog = 'REPAIR', description = 'Repairs an antibody using software tool for repairing structures (default: FoldX RepairPDB).')
    parser.add_argument('--backend', type = str, required = False, dest = 'backend', default = 'foldx', help = 'Software program to use as backend for repair (i.e. FoldX, Rosetta)', metavar = 'PKG')
    parser.add_argument('--filenames', nargs = '+', type = str, dest = 'filenames', required = True, help = 'Filenames of molecular structures to repair', metavar = 'FILENS')
    parser.add_argument('--toolpath', type = str, required=True, dest='path', help = 'Pathname of tool used for repair', metavar='PATH')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser


def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/repair/')
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

def run_repair(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)

    #### IMPLEMENTATION
    
    # > `input: path pathname of tool used for repair ` ===> namespace.path
    # `input: foldx use foldx for mutational scanning ` ===> namespace.backend
    # `input: filenames of molecular structures to repair, comma delimited` ===> namespace.filenames (this is an array of strings which are each filename)
      # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/repair/) and figs in namespace.output_dir/output/figs/
