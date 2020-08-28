import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, signal

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

def plot_predictors(filepath, output_filepath, filename, colname):
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
    if np.all(df[colname].values > 0) and np.all(df[colname+'_predicted'].values > 0):
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
    if kde_cutoff:
        plt.title('Cutoffs: positive = '+str(pos)+', negative = '+str(neg))
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
    if kde_cutoff:
        plt.title('Cutoff = '+str(cutoff))
    plt.ylabel('importances_'+colname)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_importances_plot.png', dpi=300)
