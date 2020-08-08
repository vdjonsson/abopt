import pandas as pd
import numpy as np
from io_utils import read_file, output_sparse_matrix, create_importances_dataframe, output_results, plot_predictors, plot_importances
from seqparser import seqparser
from sklearn.ensemble import RandomForestRegressor

# variables to update:
# filepath
# filename
# y_colname

filepath = '../data/'
filename = 'single_mut_effects_cleaned'

df = read_file(filepath, filename)
sparse_matrix = seqparser(df)
output_sparse_matrix(filepath, filename, sparse_matrix)

y_colname = 'bind_avg'
# offset = 331-12
y = df[y_colname].values

rf_regressor = RandomForestRegressor(random_state=0)
rf_regressor.fit(sparse_matrix, y)

predictors = rf_regressor.predict(sparse_matrix)
feature_importances = rf_regressor.feature_importances_

df_feature_importances = create_importances_dataframe(feature_importances)
output_results(filepath, filename, y_colname, df, predictors, df_feature_importances)
plot_predictors(filepath, '../figs/', filename, y_colname, log=False)
plot_importances(filepath, '../figs/', filename, y_colname)
