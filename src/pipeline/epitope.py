import argparse
import os

'''> `input: molecular structure`

    > `input: list of specific locations that includes:  WT (wild type) amino acid, chain, PDB location`

    > `input: list of chains to scan for epitopes`

    > `output: file with list of epitopes in with chain, PDB locations`'''

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'EPITOPE', description = 'Calculate epitopes of given structure')
    parser.add_argument('--filename', type = str, required = True, dest='filename', help = 'Filename of specific locations that includes: WT amino acid, chain, PDB location', metavar = 'FILEN')
    parser.add_argument('--pdb_filename', type = str, required = True, dest = 'pdb_filename', help = 'Filename including path for molecular structure (PDB file)', metavar='PDBF')
    parser.add_argument('--chains', type=str, required=True, nargs = '+', dest='chains', help = 'List of chains to scan for epitopes', metavar = 'CH')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser

def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/epitope/')
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
    pass

def run_epitope(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)
    
    ## IMPLEMENTATION

    # > `input: molecular structure` ===> namespace.pdb_filename is the filename of the PDB molecular structure

    # > `input: list of specific locations that includes:  WT (wild type) amino acid, chain, PDB location` ===> namespace.filename

    # > `input: list of chains to scan for epitopes` ===> namespace.chains
    
    # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/epitope/) and figs in namespace.output_dir/output/figs/
