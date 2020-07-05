import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# variables to update:
# filepath (line 10)
# filename (line 11)
# colname (line 14)

filepath = '../data/'
filename = 'single_mut_effects_cleaned'#'kyratsous_neutralization_data'

colname = 'bind_avg'#'REGN10989'
df = pd.read_csv(filepath+filename+'_with_'+colname+'_predictors.csv', sep=',', header=0)

plt.plot(np.sort(df.loc[:,df.columns == colname].values.flatten()), 'o')
plt.plot(np.sort(df.loc[:,df.columns == colname+'_predicted'].values.flatten()), 'o')
# plt.plot(np.log10(np.sort(df.loc[:,df.columns == colname].values.flatten())), 'o')
# plt.plot(np.log10(np.sort(df.loc[:,df.columns == colname+'_predicted'].values.flatten())), 'o')
plt.ylabel('log('+colname+')')
plt.legend([colname, 'predicted '+colname])
plt.savefig('../figs/'+filename+'_lasso_plot_'+colname+'.png')
plt.close()
