import estimator as est 
import argparse

'''
> `input: filename str of path and filename of estimator`

> `input: antibody str of antibody name`

> `input: cutoff array of [cutoffmin, cutoffmax] minimum and maximum coefficient value`

> `input optional: chain str of chain name of the molecular structure to constrain`

> `input optional: top int of number of locations to consider`

> `output: file with list of constrained features based on antibody specfic PDB locations`
'''

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'CONSTRAIN', description = 'Constrains the features of an estimator by enforcing cutoff for coefficients, using top - k features, or constraining to certain chain.')
    parser.add_argument('--filename', type=str, required=True, dest = 'filename', help = 'Filepath/filename of estimator', metavar='FILEN')
    parser.add_argument('--antibody', type = str, required=True, dest = 'antibody', help = 'Name of antibody', metavar='AB')
    parser.add_argument('--cutoff', type=float, nargs=2, dest = 'cutoff_arr', default = [None, None], help = 'Cutoff array of [cutoffmin, cutoffmax] minimum and maximum coefficient value, please specify in order: MIN MAX', metavar = 'CUTOFF')
    parser.add_argument('--chain', type=str, default=None, dest = 'chain', help = 'Chain name of molecular structure to constrain', metavar = 'CH')
    parser.add_argument('--topk', type=int, default = None, dest='top', help = 'Top number of locations to consider', metavar = 'K')
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
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')


''' Return top estimator features given an antibody name, estimator file and some cutoff '''
def top_estimator_features(antibody, estimatorfile, cutoff ):

    return

def check_args(namespace):
    # TODO
    pass

def run_constrain(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)
    
    ###### IMPLEMENTATION #######

    # > `input: filename str of path and filename of estimator` ===> namespace.filename

    # > `input: antibody str of antibody name` ===> namespace.antibody

    # > `input: cutoff array of [cutoffmin, cutoffmax] minimum and maximum coefficient value` ===> namespace.cutoff_arr (this is an array where the first element is the minimum (float value) and the second element is the maximum (float value); the default is [None, None]

    # > `input optional: chain str of chain name of the molecular structure to constrain` ===> namespace.chain (str value whose default, if not set by user, is None)

    # > `input optional: top int of number of locations to consider` ===> namespace.top (int value whose default, if not set by user, is None)
    
      # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/constrain/) and figs in namespace.output_dir/output/figs/
