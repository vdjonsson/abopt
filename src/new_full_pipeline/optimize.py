import numpy as np
import cvxpy as cp
import pandas as pd
import warnings
import argparse

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'OPTIMIZE', description = 'Uses convex combinatorial optimization to compute optimal antibody cocktail for provided virus mutants.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of fitness matrix dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of fitness matrix dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--p', type = float, required = True, dest = 'p_unmanaged', help = 'Maximum proportion of viruses not covered by the antibody cocktail', metavar = 'P_UNMANAGED')
    parser.add_argument('--l', type = float, required = True, dest = 'lmbda', help = 'Lambda value to use for penalization of number of antibodies chosen', metavar = 'LMBD')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser
    
def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/optimize/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')

def compute_antibodies(fitness_matrix, k, lmbda):
    m,n = fitness_matrix.shape
    num_unmanaged = int(m*k)
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix).clip(min=0)
    constraints = [cp.sum(c) >= 1, cp.sum_smallest((incidence_matrix@c), num_unmanaged+1) >= 1]
    objective = cp.Minimize(lmbda*cp.norm1(c)-cp.matmul(cp.sum(fitness_matrix, axis=0), c))
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value

def run_optimize(args):
    parser = setup_parser()
    namespace = parser.parse_args()
    create_output_structure(namespace.output_dir)
    
    df = pd.read_csv(namespace.filepath+namespace.filename+'.csv', sep=',', header=0, index_col=0)
    fitness_matrix = df.values
    results = compute_antibodies(fitness_matrix, namespace.p_unmanaged, namespace.lmbda)
    
    if results is None:
        warnings.warn('Maximum proportion of not managed virus is too low (no combination of antibodies can satisfy this constraint) so optimization not solved.')
    
    else:
        chosen_antibodies = df.columns[results.astype(bool)]
        with open(namespace.output_dir+'output/optimize/optimal_antibodies.txt', 'w') as f:
            f.write(','.join(chosen_antibodies))
