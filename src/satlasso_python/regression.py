import pandas as pd
import numpy as np
import satlasso_functions_cvx
import sys
from sklearn import linear_model

# variables to update:
# filepath (line 12)
# filename (line 13)
# y_colname (line 18)
# saturated (line 21)
# df_coef adjust_positions args

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

def adjust_positions(a, offset):
    new_elements=[]
    for element in a:
        pos = int(element[1:len(element)])+offset
        new_elements.append(element[0]+str(pos))
    return new_elements

filepath = '../data/'
filename = 'nussenzweig_antibody_data_cleaned_with_alignments'#'kyratsous_neutralization_data''single_mut_effects_cleaned' 'combined_single_mut_effects_cleaned_foldx_6vw1'

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)
sparse_matrix = np.genfromtxt(filepath+filename+'_sparse_matrix.csv', delimiter=',')

y_colname = sys.argv[1]
y = df[y_colname].values

saturated = True
if saturated:
    k=5
    lmbdas1_to_test = np.linspace(start=1, stop=10, num=3)
    lmbdas2_to_test = np.linspace(start=1, stop=10, num=3)
    lmbdas3_to_test = [10]# np.linspace(start=1, stop=10, num=3)
    coefficients = satlasso_functions_cvx.satlasso_solve(sparse_matrix, np.log(y).tolist(), lmbdas1_to_test, lmbdas2_to_test, lmbdas3_to_test, k)
    # predictors = np.add(np.matmul(sparse_matrix, coefficients[0:len(coefficients)-1]),[coefficients[len(coefficients)-1]]*len(sparse_matrix))
    log_predictors = np.add(np.matmul(sparse_matrix, coefficients[0:len(coefficients)-1]),[coefficients[len(coefficients)-1]]*len(sparse_matrix))
    predictors = list(map(lambda x: np.exp(x), log_predictors))
    
else:
    lasso_cv = linear_model.Lars()
    lasso_cv.fit(sparse_matrix,y)
    predictors = lasso_cv.predict(sparse_matrix)
    coefficients = lasso_cv.coef_.tolist()+[lasso_cv.intercept_]

df[y_colname+'_predicted'] = predictors
df.to_csv(filepath+filename+'_with_'+y_colname+'_predictors.csv', sep=',', index=False, header=True)

aa_positions = np.array([[s+str(i+1) for s in aalist] for i in range(0,int(len(sparse_matrix[0])/len(aalist)))]).flatten()
df_coef = pd.DataFrame(data=coefficients[0:len(coefficients)-1], index = adjust_positions(aa_positions, 0)).T
df_coef.to_csv(filepath+filename+'_'+y_colname+'_coefficients'+'.csv', sep=',', index=False, header=True)

