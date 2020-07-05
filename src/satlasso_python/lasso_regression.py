import pandas as pd
import numpy as np
import satlasso_functions
from sklearn import linear_model

# variables to update:
# filepath (line 12)
# neut_filename (line 13)
# y_colname (line 18)
# saturated (line 21)

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

filepath = '../data/'
neut_filename = 'single_mut_effects_cleaned'#'nussenzweig_antibody_data_cleaned_with_alignments''kyratsous_neutralization_data'

df = pd.read_csv(filepath+neut_filename+'.csv', sep=',', header=0)
sparse_matrix = np.genfromtxt(filepath+neut_filename+'_sparse_matrix.csv', delimiter=',')

y_colname = 'bind_avg'#'REGN10989''sars_cov_2_ic50_ngml'
y = df[y_colname].values

saturated = False
if saturated:
    k=5
    lmbdas1_to_test = [0.1, 1, 50]
    lmbdas2_to_test = [0.1]
    lmbdas3_to_test = [0.1]
    ceofficients = satlasso_functions.satlasso_solve(sparse_matrix, y, lmbdas1_to_test,lmbdas2_to_test, lmbdas3_to_test, k)
    predictors = np.add(np.matmul(sparse_matrix, coefficients[0:len(coefficients)-1]),[coefficients[len(coefficients)]]*len(sparse_matrix))
    df[y_colname+'_predicted'] = predictors
    df.to_csv(filepath+neut_filename+'_with_'+y_colname+'_predictors.csv', sep=',', index=False, header=True)
    
    aa_positions = np.array([[s+str(i+1) for s in aalist] for i in range(0,int(len(sparse_matrix[0])/len(aalist)))]).flatten()
    df_coef = pd.DataFrame(data=coefficients[0:len(coefficients)-1], index = aa_positions).T
    df_coef.to_csv(filepath+neut_filename+'_coefficients_'+y_colname+'.csv', sep=',', index=False, header=True)
else:
    lasso_cv = linear_model.LassoCV(cv=5)
    lasso_cv.fit(sparse_matrix,y)
    predictors = lasso_cv.predict(sparse_matrix)
    coefficients = lasso_cv.coef_.tolist()+[lasso_cv.intercept_]
    df[y_colname+'_predicted'] = predictors
    df.to_csv(filepath+neut_filename+'_with_'+y_colname+'_predictors.csv', sep=',', index=False, header=True)
    
    aa_positions = np.array([[s+str(i+1) for s in aalist] for i in range(0,int(len(sparse_matrix[0])/len(aalist)))]).flatten()
    df_coef = pd.DataFrame(data=coefficients[0:len(coefficients)-1], index = aa_positions).T
    df_coef.to_csv(filepath+neut_filename+'_coefficients_'+y_colname+'.csv', sep=',', index=False, header=True)
