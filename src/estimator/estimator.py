import pandas as pd
import numpy as np
import sys
import os
import argparse
from SatLasso import SatLasso, SatLassoCV
from io_utils import read_file, output_sparse_matrix, create_coefs_dataframe, output_results, output_mapped_coefs, output_opt_lambdas
from plot_utils import plot_predictors, plot_coefs
from seqparser import seqparser, map_coefs

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'ESTIMATOR', description = 'Parse amino acid sequence and run SatLasso for variable selection.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--s', required = False, default = 'max', dest = 'saturation', help = 'Saturation value to use for SatLasso(CV): can be float or {"max", "mode"}', metavar = 'SAT')
    parser.add_argument('--m', action = 'store_const', const = True, required = False, dest = 'map_back', default = False, help = 'Use flag to map coefficients back to individual amino acid sequences')
    parser.add_argument('--y', type = str, required = True, dest = 'y_colname', help = 'Name of column for y values', metavar = 'Y')
    parser.add_argument('--cv', type = int, required = False, dest = 'cv', default = 0, help = 'Cross-validation value: use int for SatLassoCV with specified number of folds; otherwise SatLasso (no CV) used', metavar = 'CV')
    parser.add_argument('--t', choices = ['log10', 'ln', 'norm'], required = False, dest = 'transform', metavar = 'TRNSFM', default = None, help = 'Transform to use on y values to feed into cross-validation: can be {"ln", "log10", "norm"}')
    parser.add_argument('--l1', type = float, required = True, dest = 'lambda1', help = 'Lambda 1 value; or if using CV, start value for lambda 1', metavar = 'LMBD1')
    parser.add_argument('--l2', type = float, required = True, dest = 'lambda2', help = 'Lambda 2 value; or if using CV, start value for lambda 2', metavar = 'LMBD2')
    parser.add_argument('--l3', type = float, required = True, dest = 'lambda3', help = 'Lambda 3 value; or if using CV, start value for lambda 3', metavar = 'LMBD3')
    parser.add_argument('--n', type = int, required = False, default = 10, dest = 'n_lambdas', help = 'Number of lambdas to use for grid search in CV; ignored if not using CV', metavar = 'N_LMBDS')
    parser.add_argument('--r', type = float, required = False, default = 10, dest = 'range', help = 'Range to use for grid search of lambdas in CV; ignored if not using CV', metavar = 'RANGE')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser
    
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

def check_dataframe(df, map_back):
    if 'sequence' not in df.columns:
        raise ValueError('Dataframe object must have column with sequences named "sequence".')
    if np.unique(map(len, df['sequence'].values)).size > 1:
        raise ValueError('Sequence column must have sequences of equal sizes.')
    if 'heavy_chain_aligned' not in df.columns and map_back:
        raise ValueError('Dataframe object must have column with aligned heavy chain sequences named "heavy_chain_aligned" if map back coefficients set to True.')
    if 'light_chain_aligned' not in df.columns and map_back:
        raise ValueError('Dataframe object must have column with aligned light chain sequences named "light_chain_aligned" if map back coefficients set to True.')
    if 'antibody_id' not in df.columns and map_back:
        raise ValueError('Dataframe object must have column identifying each sequence named "antibody_id" if map back coefficients set to True.')
    
def main():
    parser = setup_parser()
    namespace = parser.parse_args()
    create_output_structure(namespace.output_dir)
    
    df = read_file(namespace.filepath, namespace.filename)
    check_dataframe(df, namespace.map_back)
    sparse_matrix = seqparser(df, 'sequence')
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
        mapped_coefs = map_coefs(df, df_coefs, 'heavy_chain_aligned', 'light_chain_aligned', 'antibody_id', 'sequence')
        output_mapped_coefs(namespace.output_dir+'output/estimator/data/', namespace.filename, mapped_coefs)
    
    plot_predictors(namespace.output_dir+'output/estimator/data/', namespace.output_dir+'output/estimator/figs/', namespace.filename, namespace.y_colname)
    plot_coefs(namespace.output_dir+'output/estimator/data/', namespace.output_dir+'output/estimator/figs/', namespace.filename, namespace.y_colname)

if __name__ == "__main__":
    main()
