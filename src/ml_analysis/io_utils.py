import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seqparser import adjust_positions
from scipy import stats, signal

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

def read_file(filepath, filename):
    """
    Reads pandas DataFrame from CSV file.
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
        
    Returns
    ----------
    df : pandas DataFrame
        DataFrame from file
    """
    
    df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)
    return df

def output_sparse_matrix(filepath, filename, sparse_matrix):
    """
    Output numpy array to CSV file.
    
    Note: Intended to be used as utility function for I/O of sparse matrix generated /
        by seqparser.
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
    sparse_matrix : numpy array
        Array to save in file
    """
    
    np.savetxt(filepath+filename+'_sparse_matrix.csv', sparse_matrix, delimiter=',')

def read_sparse_matrix(filepath, filename):
    """
    Reads/generates numpy array from CSV file.
    
    Note: Intended to be used as utility function for I/O of sparse matrix generated /
        by seqparser.
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
        
    Returns
    ----------
    sparse_matrix : numpy array
        Numpy array from file
    """
    
    sparse_matrix = np.genfromtxt(filepath+filename+'_sparse_matrix.csv', delimiter=',')
    return sparse_matrix

def create_coefs_dataframe(coefs, offset=0):
    """
    Create dataframe from coefficients with index : amino acid positions.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    coefs : ndarray of shape (n_coefficients,)
        Coefficient values
    offset : int, default = 0
        Integer change to pass to adjust_positions (see function adjust_positions in seqparser)
        
    Returns
    ----------
    df : pandas DataFrame
        Coefficient dataframe
    """
    
    aa_positions = np.array([[s+str(i) for s in aalist] for i in range(0,int(len(coefs)/len(aalist)))]).flatten()
    df = pd.DataFrame(data=coefs, index = adjust_positions(aa_positions, offset), columns = ['coefficients'])
    return df

def create_importances_dataframe(importances, offset=0):
    """
    Create dataframe from importances with index : amino acid positions.
    
    Note: Intended to be used with feature importances output from /
        RandomForestRegressor
    
    Parameters
    ----------
    importances : ndarray of shape (n_importances,)
        Feature importance values
    offset : int, default = 0
        Integer change to pass to adjust_positions (see function adjust_positions in seqparser)
        
    Returns
    ----------
    df : pandas DataFrame
        Importance dataframe
    """
    aa_positions = np.array([[s+str(i) for s in aalist] for i in range(0,int(len(importances)/len(aalist)))]).flatten()
    df = pd.DataFrame(data=importances, index = adjust_positions(aa_positions, offset), columns = ['importances'])
    return df

    
def output_results(filepath, filename, colname, df, predictors, df_coefs):
    """
    Utility function to output results of regression to CSV files.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
    colname : str
        Name of column that results predict
    df : pandas DataFrame
        Metadata
    predictors : ndarray of shape (n_predictors,)
        Predicted result values
    df_coefs : pandas DataFrame
        Coefficient dataframe
    """
    
    df[colname+'_predicted'] = predictors
    df.to_csv(filepath+filename+'_with_predictors.csv', sep=',', header=True, index=False)
    df_coefs.to_csv(filepath+filename+'_coefficients.csv', sep=',', header=True, index=True)

def plot_predictors(filepath, output_filepath, filename, colname, log=True):
    """
    Plot and output predicted values and true values.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    filepath : str
        Name of path
    output_filepath : str
        Name of path to output path
    filename : str
        Name of file
    colname : str
        Name of column that results predict
    log : bool, default = True
        Whether to plot on log scale
    """
    
    df = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
    df.sort_values(by=[colname], inplace=True)
    if log:
        plt.plot(np.log10(df.loc[:,df.columns == colname].values.flatten()), 'o')
        plt.plot(np.log10(df.loc[:,df.columns == colname+'_predicted'].values.flatten()), 'o')
        plt.ylabel('log('+colname+')')
        plt.legend(['log '+colname, 'log predicted '+colname])
    else:
        plt.plot(df.loc[:,df.columns == colname].values.flatten(), 'o')
        plt.plot(df.loc[:,df.columns == colname+'_predicted'].values.flatten(), 'o')
        plt.ylabel(colname)
        plt.legend([colname, 'predicted '+colname])
    plt.savefig(output_filepath+filename+'_predictors_plot.png', dpi=300)
    plt.close()

def coefficient_cutoff(non_zero_coeff):
    """
    Utility function for adjusting coefficients returned by SatLasso /
        in order to remove coefficients with very low values. Calculates /
        thresholds for cutoff based on Gaussian kernel density estimation /
        of coefficient values.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    non_zero_coeff: ndarray of shape (n_coefficients,)
        Coefficient values
        
    Returns
    ----------
    negative_cutoff : float
        Threshold for negative coefficients
    positive_cutoff : float
        Threshold for positive coefficients
    """
    
    kde = stats.gaussian_kde(non_zero_coeff)
    x = np.linspace(non_zero_coeff.min(), non_zero_coeff.max(), num = len(np.unique(non_zero_coeff)))
    y = kde.evaluate(x)
    valleys = x[signal.argrelextrema(y, np.less)]
    negative_cutoff = max([n for n in valleys if n<0])
    positive_cutoff = min([n for n in valleys if n>0])
    return (negative_cutoff, positive_cutoff)

def plot_coefs(filepath, output_filepath, filename, colname, kde_cutoff = True):
    """
    Plot and output bar plot of coefficients.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    filepath : str
        Name of path
    output_filepath : str
        Name of path to output path
    filename : str
        Name of file
    colname : str
        Name of column that results predict
    kde_cutoff : bool, default = True
        Whether to use Gaussian kernel density estimation to calculate threshold /
            for coefficients
    """
    
    df = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
    non_zero_coefs = df.loc[df.coefficients != 0]
    if kde_cutoff:
        neg, pos = coefficient_cutoff(non_zero_coefs.coefficients.values)
        non_zero_coefs = non_zero_coefs.loc[np.logical_or(non_zero_coefs.coefficients.values < neg, non_zero_coefs.coefficients.values > pos)]
    plt.figure(figsize=(20,5), dpi=300)
    plt.bar(range(0, len(non_zero_coefs)), non_zero_coefs.coefficients.values, tick_label= non_zero_coefs.index)
    plt.xticks(rotation=90)
    plt.ylabel('coeff_'+colname)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_coefficients_plot.png', dpi=300)
    
def importance_cutoff(non_zero_importances):
    """
    Utility function for adjusting feature importances returned by RandomForestRegressor /
        in order to remove feature importances with very low values. Calculates /
        thresholds for cutoff based on Gaussian kernel density estimation /
        of importance values.
    
    Note: Intended to be used with feature importance output from RandomForestRegressor
    
    Parameters
    ----------
    non_zero_importances: ndarray of shape (n_importances,)
        Feature importance values
        
    Returns
    ----------
    cutoff : float
        Threshold for feature importances
    """
    
    kde = stats.gaussian_kde(non_zero_importances)
    x = np.linspace(non_zero_importances.min(), non_zero_importances.max(), num = len(np.unique(non_zero_importances)))
    y = kde.evaluate(x)
    valleys = x[signal.argrelextrema(y, np.less)]
    cutoff = valleys[0]
    return cutoff

def plot_importances(filepath, output_filepath, filename, colname, kde_cutoff = True):
    """
    Plot and output bar plot of feature importances.
    
    Note: Intended to be used with importance output from RandomForestRegressor
    
    Parameters
    ----------
    filepath : str
        Name of path
    output_filepath : str
        Name of path to output path
    filename : str
        Name of file
    colname : str
        Name of column that results predict
    kde_cutoff : bool, default = True
        Whether to use Gaussian kernel density estimation to calculate threshold /
            for feature importances
    """
    
    df = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
    non_zero_coefs = df.loc[df.importances != 0]
    if kde_cutoff:
        cutoff = importance_cutoff(non_zero_coefs.importances.values)
        non_zero_coefs = non_zero_coefs.loc[non_zero_coefs.importances.values > cutoff]
    plt.figure(figsize=(20,5), dpi=300)
    plt.bar(range(0, len(non_zero_coefs)), non_zero_coefs.importances.values, tick_label= non_zero_coefs.index)
    plt.xticks(rotation=90)
    plt.ylabel('importances_'+colname)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_importances_plot.png', dpi=300)

def output_mapped_coefs(filepath, filename, coefs_df):
    """
    Output coefficients mapped back to original positions in sequences.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
    coefs_df : pandas DataFrame
        Coefficients mapped back to amino acid in each heavy/light chain sequence with /
            MultiIndex : identifying_name, location, chain, amino_acid and columns: /
            wild_type and coefficient
    """
    
    coefs_df.to_csv(filepath+filename+'_mapped_coefficients.csv', sep=',', header=True, index=True)
