import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

def create_coefs_dataframe(coefs):
    """
    Create dataframe from coefficients with index : amino acid positions.
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV).
    
    Parameters
    ----------
    coefs : ndarray of shape (n_coefficients,)
        Coefficient values
        
    Returns
    ----------
    df : pandas DataFrame
        Coefficient dataframe
    """
    
    aa_positions = np.array([[s+str(i) for s in aalist] for i in range(0,int(len(coefs)/len(aalist)))]).flatten()
    df = pd.DataFrame(data=coefs, index = aa_positions, columns = ['coefficients'])
    return df

def create_importances_dataframe(importances):
    """
    Create dataframe from importances with index : amino acid positions.
    
    Note: Intended to be used with feature importances output from /
        RandomForestRegressor
    
    Parameters
    ----------
    importances : ndarray of shape (n_importances,)
        Feature importance values
        
    Returns
    ----------
    df : pandas DataFrame
        Importance dataframe
    """
    aa_positions = np.array([[s+str(i) for s in aalist] for i in range(0,int(len(importances)/len(aalist)))]).flatten()
    df = pd.DataFrame(data=importances, index = aa_positions, columns = ['importances'])
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

def output_opt_lambdas(filepath, filename, lambda1, lambda2, lambda3):
    """
    Output optimal lambda values
    
    Parameters
    ----------
    filepath : str
        Name of path
    filename : str
        Name of file
    lambda1 : float
        optimal lambda 1 value
    lambda2 : float
        optimal lambda 2 value
    lambda3 : float
        optimal lambda 3 value
    """
    
    df = pd.DataFrame(data=[lambda1, lambda2, lambda3], index=['lambda1', 'lambda2', 'lambda3'], columns=['optimal'])
    df.to_csv(filepath+filename+'_optimal_lambdas.csv', sep=',', index=True, header=True)
