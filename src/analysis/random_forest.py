import pandas as pd
import numpy as np
from io_utils import read_file, output_sparse_matrix, create_importances_dataframe, output_results
from seqparser import seqparser
from plot_utils import plot_predictors, plot_importances
from sklearn.ensemble import RandomForestRegressor

# variables to update:
# filepath
# filename
# y_colname

def adjust_positions(a, offset=0):
    """
    Utility function to adjust position notation in list of amino acid positions with form /
        <amino acid><position>, for example: A34.

    Parameters
    ----------
    a : ndarray of shape (n_positions,)
        Amino acid positions to adjust
    offset : int, default = 0
        Integer change in each amino acid position

    Returns
    ----------
    new_elements : ndarray of shape (n_positions,)
        Amino acid positions adjusted
    """

    if offset == 0:
        return a
        
    new_elements=[]
    for element in a:
        pos = int(element[1:len(element)])+offset
        new_elements.append(element[0]+str(pos))
    return new_elements

filepath = '../data/'
output_filepath = '../output/estimator/data/'
output_fig_filepath = '../output/estimator/figs/'
filename = 'single_mut_effects_cleaned'
seq_col = 'sequence'

df = read_file(filepath, filename)
sparse_matrix = seqparser(df, seq_col)
output_sparse_matrix(output_filepath, filename, sparse_matrix)

y_colname = 'bind_avg'
offset = 331-12
y = df[y_colname].values

rf_regressor = RandomForestRegressor(random_state=0)
rf_regressor.fit(sparse_matrix, y)

predictors = rf_regressor.predict(sparse_matrix)
feature_importances = rf_regressor.feature_importances_

df_feature_importances = create_importances_dataframe(feature_importances)
df_feature_importances.index = adjust_positions(df_feature_importances.index.values, offset)
output_results(output_filepath, filename, y_colname, df, predictors, df_feature_importances)
plot_predictors(output_filepath, output_fig_filepath, filename, y_colname)
plot_importances(output_filepath, output_fig_filepath, filename, y_colname)
