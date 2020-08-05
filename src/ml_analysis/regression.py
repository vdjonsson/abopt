import pandas as pd
import numpy as np
from SatLasso import SatLasso, SatLassoCV
from io_utils import read_file, output_sparse_matrix, create_coefs_dataframe, output_results, plot_predictors, plot_coefs, output_mapped_coefs
from seqparser import seqparser, align_coefs

# variables to update:
# filepath
# filename
# y_colname
# lambdas to test
# offset

filepath = '../data/'
filename = 'NeutSeqData_VH3-53_66_aligned'#'NeutSeqData_C002-215_cleaned_aligned'
heavy_chain = 'VH or VHH'#'igh_vdj_aa'
light_chain = 'VL'#'igl_vj_aa'
id_col = 'Name'#'antibody_id'
y_colname = 'IC50_ngml'#'sars_cov_2_ic50_ngml'
offset = 0

df = read_file(filepath, filename)
sparse_matrix = seqparser(df)
output_sparse_matrix(filepath, filename, sparse_matrix)

y = df[y_colname].values.astype(float)

lmbdas1_to_test = np.linspace(start = 1, stop = 10, num=5)
lmbdas2_to_test = np.linspace(start = 1, stop = 10, num=5)
lmbdas3_to_test = np.linspace(start = 5, stop = 10, num=5)
satlassoCV = SatLassoCV(lambda_1s = lmbdas1_to_test, lambda_2s = lmbdas2_to_test, lambda_3s = lmbdas3_to_test, saturation='mode', cv=3)
satlassoCV.fit(sparse_matrix, y)
coefficients = satlassoCV.coef_
predictors = satlassoCV.predict(sparse_matrix)
print(satlassoCV.lambda_1_, satlassoCV.lambda_2_, satlassoCV.lambda_3_)
# log_predictors = satlassoCV.predict(sparse_matrix)
# predictors = list(map(lambda x: np.exp(x), log_predictors))

df_coefs = create_coefs_dataframe(coefficients, offset)
output_results(filepath, filename, y_colname, df, predictors, df_coefs)
mapped_coefs = map_coefs(df, df_coefs, heavy_chain, light_chain, id_col)
output_mapped_coefs(filepath, filename, mapped_coefs)
plot_predictors(filepath, '../figs/', filename, y_colname)
plot_coefs(filepath, '../figs/', filename, y_colname)
