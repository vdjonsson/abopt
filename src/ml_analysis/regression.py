import pandas as pd
import numpy as np
import sys
import os
import warnings
import argparse
from SatLasso import SatLasso, SatLassoCV
from io_utils import read_file, output_sparse_matrix, create_coefs_dataframe, output_results, output_mapped_coefs, output_opt_lambdas
from plot_utils import plot_predictors, plot_coefs
from seqparser import seqparser, map_coefs

def setup_parser():
    parser = argparse.ArgumentParser(prog = '<PROGRAM NAME>', description = 'Parse amino acid sequence and run SatLasso.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'filepath of data', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'filename of data', metavar = 'FILEN')
    parser.add_argument('--seq', type = str, required = True, dest = 'seq_col', help = 'name of column in dataframe with sequences', metavar = 'SEQCOL')
    parser.add_argument('--sat', required = False, default = 'max', dest = 'saturation', help = 'saturation value to use for SatLasso(CV); can be float or {"max", "mode"}', metavar = 'SAT')
    parser.add_argument('--type', choices = ['virus', 'antibody'], required = False, default = 'antibody', dest = 'type', help = 'type of structure characterized by AA sequences (virus or antibody)')
    parser.add_argument('--hc', type = str, required = False, dest = 'heavy_chain_col', help = 'name of aligned heavy chain column in dataframe', metavar = 'HEAVY', default = None)
    parser.add_argument('--lc', type = str, required = False, dest = 'light_chain_col', help = 'name of aligned light chain column in dataframe', metavar = 'LIGHT', default = None)
    parser.add_argument('--a', action = 'store_const', const = True, required = False, dest = 'map_back', default = False, help = 'whether sequences are aligned and need to be mapped back')
    parser.add_argument('--id', type = str, required = True, dest = 'id_col', help = 'name of identifying column for antibody in dataframe', metavar = 'ID')
    parser.add_argument('--y', type = str, required = True, dest = 'y_colname', help = 'name of column for y values', metavar = 'Y')
    parser.add_argument('--cv', type = int, required = False, dest = 'cv', default = 0, help = 'whether to use cross validation, use int for SatLassoCV with specified number of folds, otherwise SatLasso (no CV) used', metavar = 'CV')
    parser.add_argument('--t', choices = ['log10', 'ln', 'norm'], required = False, dest = 'transform', metavar = 'TRNSFM', default = None, help = 'whether to transform data and how to transform data')
    parser.add_argument('--l1', type = float, required = True, dest = 'lambda1', help = 'lambda 1 value, or if using CV, start value for lambda 1', metavar = 'LMBD1')
    parser.add_argument('--l2', type = float, required = True, dest = 'lambda2', help = 'lambda 2 value, or if using CV, start value for lambda 2', metavar = 'LMBD2')
    parser.add_argument('--l3', type = float, required = True, dest = 'lambda3', help = 'lambda 3 value, or if using CV, start value for lambda 3', metavar = 'LMBD3')
    parser.add_argument('--n', type = int, required = False, default = 10, dest = 'n_lambdas', help = 'number of lambdas to use for grid search in CV; ignored if not using CV', metavar = 'N_LMBDS')
    parser.add_argument('--r', type = float, required = False, default = 10, dest = 'range', help = 'range to use for grid search in CV; ignored if not using CV', metavar = 'RANGE')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'output directory', metavar = 'OUT')
    parser.add_argument('--v', action = 'store_const', const = True, required = False, default = False, dest = 'verbose', help = 'verbose')
    return parser

def check_arguments(namespace):
    if namespace.map_back and namespace.type == 'antibody' and namespace.heavy_chain_col is None and namespace.light_chain_col is None:
        warnings.warn('If mapping back to aligned antibodies and no heavy chain or light chain column is provided, the sequence column will be used to map back coefficients (no chain specific information encoded).')
    if namespace.map_back and namespace.type == 'virus' and namespace.heavy_chain_col is not None and namespace.light_chain_col is not None:
        warnings.warn('If mapping back to aligned viruses, heavy chain and light chain arguments will be ignored and sequence column will be used.')
    return True
    
def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/estimator/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/estimator/data')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/estimator/figs')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')
    
def main():
    parser = setup_parser()
    namespace = parser.parse_args()
    check_arguments(namespace)
    create_output_structure(namespace.output_dir)
    
    df = read_file(namespace.filepath, namespace.filename)
    sparse_matrix = seqparser(df, namespace.seq_col)
    output_sparse_matrix(namespace.output_dir+'output/estimator/data/', namespace.filename, sparse_matrix)
    y = df[namespace.y_colname].values.astype(float)
    
    if namespace.transform == 'ln':
        y = np.log(y)
    elif namespace.transform == 'log10':
        y = np.log10(y)
    
    if not namespace.cv:
        satlasso = SatLasso(lambda_1 = namespace.lambda1, lambda_2 = namespace.lambda2, lambda_3 = namespace.lambda3, saturation = namespace.saturation, normalize = (namespace.transform == 'norm'))
            
    else:
        lambda_1s_grid = np.linspace(start = namespace.lambda1, stop = namespace.lambda1+namespace.range, num = namespace.n_lambdas)
        lambda_2s_grid = np.linspace(start = namespace.lambda2, stop = namespace.lambda2+namespace.range, num = namespace.n_lambdas)
        lambda_3s_grid = np.linspace(start = namespace.lambda3, stop = namespace.lambda3+namespace.range, num = namespace.n_lambdas)
        if isinstance(namespace.cv, bool):
            satlasso = SatLassoCV(lambda_1s = lambda_1s_grid, lambda_2s = lambda_2s_grid, lambda_3s = lambda_3s_grid, saturation = namespace.saturation, normalize = (namespace.transform == 'norm'))
        else:
            satlasso = SatLassoCV(lambda_1s = lambda_1s_grid, lambda_2s = lambda_2s_grid, lambda_3s = lambda_3s_grid, saturation = namespace.saturation, cv = namespace.cv)
    
    satlasso.fit(sparse_matrix, y)
    coefficients = satlasso.coef_
    
    if namespace.transform == 'ln':
        log_predictors = satlasso.predict(sparse_matrix)
        predictors = list(map(lambda x: np.exp(x), log_predictors))
    elif namespace.transform == 'log10':
        log_predictors = satlasso.predict(sparse_matrix)
        predictors = list(map(lambda x: 10**x, log_predictors))
    else:
        predictors = satlasso.predict(sparse_matrix)
    
    df_coefs = create_coefs_dataframe(coefficients)
    output_results(namespace.output_dir+'output/estimator/data/', namespace.filename, namespace.y_colname, df, predictors, df_coefs)
    
    if namespace.cv:
        output_opt_lambdas(namespace.output_dir+'output/estimator/data/', namespace.filename, satlasso.lambda_1_, satlasso.lambda_2_, satlasso.lambda_3_)
    
    if namespace.map_back:
        mapped_coefs = map_coefs(df, df_coefs, namespace.heavy_chain_col, namespace.light_chain_col, namespace.id_col, namespace.seq_col)
        output_mapped_coefs(namespace.output_dir+'output/estimator/data/', namespace.filename, mapped_coefs)
    
    plot_predictors(namespace.output_dir+'output/estimator/data/', namespace.output_dir+'output/estimator/figs/', namespace.filename, namespace.y_colname)
    plot_coefs(namespace.output_dir+'output/estimator/data/', namespace.output_dir+'output/estimator/figs/', namespace.filename, namespace.y_colname)

if __name__ == "__main__":
    main()
