import pandas as pd
import numpy as np
import seaborn as sb
import sys
from scipy import stats, signal
import logomaker
import matplotlib.pyplot as plt

# variables to update:
# method (line 14)
# filepath (line 10)
# filename (line 11)
# colname (line 14)

def coefficient_cutoff(non_zero_coeff):
    kde = stats.gaussian_kde(non_zero_coeff)
    x = np.linspace(non_zero_coeff.min(), non_zero_coeff.max(), num = len(np.unique(non_zero_coeff)))
    y = kde.evaluate(x)
    valleys = x[signal.argrelextrema(y, np.less)]
    negative_cutoff = max([n for n in valleys if n<0])
    positive_cutoff = min([n for n in valleys if n>0])
    return (negative_cutoff, positive_cutoff)

method = 'satlasso'
filepath = '../data/'
filename = 'nussenzweig_antibody_data_cleaned_with_alignments'#'single_mut_effects_cleaned''kyratsous_neutralization_data''combined_single_mut_effects_cleaned_foldx_6vw1'

colname = sys.argv[1]#'sars_cov_2_ic50_ngml'#'REGN10986''bind_avg'
df = pd.read_csv(filepath+filename+'_with_'+colname+'_predictors.csv', sep=',', header=0)
df.sort_values(by=[colname], inplace=True)

#barplot
# plt.plot(df.loc[:,df.columns == colname].values.flatten(), 'o')
# plt.plot(df.loc[:,df.columns == colname+'_predicted'].values.flatten(), 'o')
plt.plot(np.log10(df.loc[:,df.columns == colname].values.flatten()), 'o')
plt.plot(np.log10(df.loc[:,df.columns == colname+'_predicted'].values.flatten()), 'o')
plt.ylabel('log('+colname+')')
plt.legend([colname, 'predicted '+colname])
plt.savefig('../figs/'+filename+'_'+method+'_plot_'+colname+'.png')
plt.close()

df = pd.read_csv(filepath+filename+'_'+colname+'_coefficients.csv', sep=',', header=0)
non_zero_coeff = df.iloc[0,:].values[df.iloc[0,:].values!=0]

# optional cutoff
# neg, pos = coefficient_cutoff(non_zero_coeff)
# non_zero_coeff = df.iloc[0,:].values[np.logical_or(df.iloc[0,:].values < neg, df.iloc[0,:].values > pos)]

plt.figure(figsize=(20,5), dpi=300)
plt.bar(range(0, len(non_zero_coeff)), non_zero_coeff, tick_label = df.columns[df.iloc[0,:].values!=0])
plt.xticks(rotation=90)
plt.ylabel(method+'_coeff_'+colname)
plt.tight_layout()
plt.savefig('../figs/'+filename+'_bar_'+method+'_coeff_'+colname+'.png')
plt.close()

# # heatmap
# aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

# aa_heatmap = pd.DataFrame(index=aalist)

# for i in range(0, int(len(df.columns)/len(aalist))):
#     index = df.columns[i*len(aalist)][1:]
#     aa_heatmap[index] = df.iloc[0,:][i*len(aalist):(i+1)*len(aalist)].values

# logo_colors = sb.color_palette('husl', len(aalist))
# logo_heatmap = aa_heatmap.iloc[:,11:213].T
# logo_heatmap.index = logo_heatmap.index.astype(int)
# plt.figure(figsize=(15,5), dpi=300)
# logo = logomaker.Logo(logo_heatmap, shade_below=0.5, fade_below=0.5, show_spines = False, color_scheme = dict(zip(aalist, logo_colors)), figsize=(15,5))
# logo.ax.set_ylabel('coefficient')
# plt.tight_layout()
# plt.savefig('../figs/'+filename+'_logo_plot_'+colname+'.png', dpi='figure')
# plt.close()

# plt.figure(figsize=(15,5), dpi=300)
# sb.heatmap(aa_heatmap.iloc[:,11:213])
# plt.savefig('../figs/'+filename+'_coeff_heatmap_'+colname+'.png', dpi='figure')
# plt.close()

# # plot violin plot (?)
# plt.figure(figsize=(15,5), dpi=300)
# sb.violinplot(data = aa_heatmap.iloc[:,11:213])
# plt.xticks(rotation = 90, fontsize=4)
# plt.ylabel('coefficient')
# plt.tight_layout()
# plt.savefig('../figs/'+filename+'_violin_plot_'+colname+'_coefficients.png', dpi='figure')
# plt.close()

# df = pd.read_csv(filepath+filename+'_with_'+colname+'_predictors.csv', sep=',', header=0)

# aa_heatmap = pd.DataFrame(data = np.array([np.nan]*len(aalist)*len(np.arange(331,532))).reshape((len(aalist),len(np.arange(331,532)))), index=aalist, columns=np.arange(331,532))

# for i in range(0, len(df)):
#     aa_heatmap[df.site_SARS2.values[i]][df.mutant.values[i]] = df[colname+'_predicted'].values[i]

# plt.figure(figsize=(15,5), dpi=300)
# sb.heatmap(aa_heatmap, mask = aa_heatmap.isnull())
# plt.savefig('../figs/'+filename+'_heatmap_'+colname+'_predicted.png', dpi='figure')
# plt.close()
