import sys
sys.path.insert(1, '/Users/vjonsson/Google Drive/data/repository/abopt-private/src/pipeline')
 
import cocktail
import numpy as np 

import itertools as it 
import pandas as pd

''' Move this out ''' 
import seaborn as sb 
import matplotlib.pyplot as plt 
import colors 

cocktail.create_output_structure(output_dir ='../../')

coverage = np.linspace(0.06, 0.55, num = 20)
coverage  = np.round(coverage, 2) 

#coverage = [0.2, 0.3, 0.4]

ligand1_fitness = '../../output/merge/rbd_ab_fitness_opt.csv'
ligand2_fitness = '../../output/merge/rbd_ace2_fitness.csv'


''' Log step gamma1 between 2000 and 10 '''

''' Re run everything with num= 50 '''
numsims = 10
gamma1 = np.logspace(start=0 , stop=15, num=30, endpoint=True, base=2.0)
gamma2 = np.logspace(start=0 , stop=15, num=30, endpoint=True, base=2.0)

noise_sims = 10
cocktail.run_simulations(ligand1_fitness, ligand2_fitness, coverage,gamma1, gamma2, noise_sims=noise_sims)    

exit()
''' Analyze post running sims '''

figname = 'graph_scatter_gamma'
data = pd.read_csv('../../output/cocktail/cocktails_allsims.csv')


''' Uniquely name these combinations ''' 

print(data.columns)
uc = data.iloc[:,0:-4]
uc = uc.drop_duplicates()
uc = uc.astype(int).astype(str)
uc['new'] = uc.apply(''.join, axis=1)


print(uc.head())

cnames = ['C' + str(x) for x in range(len(uc))]
names_dict = dict(zip(uc.new.values,cnames)) 

data['new'] = data.iloc[:,0:-4].astype(int).astype(str).apply(''.join, axis=1)
datanames = [names_dict[c] for c in data.new.values]
data['cocktail'] = datanames 

''' Calculate infective virus coverage ''' 

abdf, vdf = cocktail.import_fitness_landscapes(ligand1_fitness, ligand2_fitness)

fitness_matrix = abdf.values
virus_matrix = vdf.values
incidence_matrix = np.sign(fitness_matrix).clip(min=0)
virus_incidence_matrix = np.sign(virus_matrix).clip(min=0)

vals, vals2 = [] , []

''' Get list of unique cocktails'''

uc = data.iloc[:,0:-7] 
uc ['cocktail'] = data.cocktail.values 
uc  = uc.drop_duplicates()

for i in range(len(uc)):
    c = uc.iloc[i,:-1].values 
    vals.append(np.matmul(np.matmul(incidence_matrix,c),virus_incidence_matrix)[0])
    vals2.append(np.matmul(np.sum(fitness_matrix, axis=0), c))

vals_dict = dict(zip(uc.cocktail.values, vals))
vals2_dict = dict(zip(uc.cocktail.values, vals2))

data['num_abs'] = data.num_abs.astype(int)
data['inf_cov'] = [vals_dict[c] for c in data.cocktail.values]
data['neut_cov'] = [vals2_dict[c] for c in data.cocktail.values]
data['opt_cov'] = data.neut_cov - data.inf_cov


''' Graph all '''
plt.figure(figsize=(3.25,2.5))
sb.set(context='paper', style='ticks') 
sb.scatterplot(data=data, x='num_abs', y='opt_cov', hue='cocktail',legend=False)
plt.xlabel('Abs')
plt.ylabel('Infective virus covered')
plt.tight_layout()
plt.savefig('../../output/figs/cocktail_abs.png', dpi=300, transparent=True)
plt.show()

labels = ['$\gamma_1$','$\gamma_2$']
cols = ['gamma1','gamma2']
dict = dict(zip(cols, labels))
data = data.rename(columns=dict) 

optvals = pd.DataFrame()
for c in data['cov'].unique(): 
    datass = data.loc[data['cov'] == c][['cov', 'cocktail','inf_cov', 'neut_cov', 'num_abs', 'opt_cov']]
    minabs = datass.num_abs.min()
    tmp = datass.loc[datass.num_abs == minabs]
    tmp = tmp.drop_duplicates()

    '''Find cocktails with minimum number of antibodies and max infective virus coverage''' 
    ''' Either inf_cov or opt_cov '''
    ''' opt_cov '''
    minval = tmp.opt_cov.min()
    tmp = tmp.loc[tmp.opt_cov == minval]

    ''' Inf_cov '''
    # maxval = tmp.inf_cov.max()
    # tmp = tmp.loc[tmp.inf_cov == maxval]
    optvals = pd.concat([optvals, tmp])

print (optvals)

optvals = optvals.set_index('cov')
pal = sb.color_palette("cubehelix", len(optvals.cocktail.values))

''' Tradeoffs num abs vs. inf virus coverage  ''' 
size = optvals.num_abs.values**1.5*7 
plt.figure(figsize=(3.25,2.5))
sb.set(context='paper', style='ticks') 
sb.scatterplot(data=optvals, x='cov', y='inf_cov', hue='cocktail', s=size,  color=pal,legend=False)
sb.lineplot(data=optvals, x='cov', y='inf_cov',color = 'grey', alpha=0.7, lw=0.3)
plt.xlabel('% Minimum coverage')
plt.ylabel('Infective virus covered')
#plt.legend(ncol=3)
plt.tight_layout()
plt.savefig('../../output/figs/cocktail_optimal.png', dpi=300, transparent=True)
plt.show()

optvals['ccov'] = optvals.index.astype(str) + optvals.cocktail
data['ccov'] = data['cov'].astype(str) + data.cocktail
isopt = []

for i in range(len(data)):
    isopt.append(data.iloc[i, :].ccov in optvals.ccov.values)

print(np.sum(isopt))
data['isopt'] = isopt


''' Graph heatmap of optimal combinations per coverage ''' 

dataopt = data.loc[data.isopt == True]
dataopt = dataopt.drop(['ccov','new','$\gamma_1$', '$\gamma_2$'],axis=1)
dataopt = dataopt.drop_duplicates('cocktail')

cocktails = dataopt.cocktail.unique()

cocktail.graph_cocktail_mixes_data(dataopt)

''' Facet Grid'''
height = 1.75
aspect = 0.75

ymax = 80
ms = 4

''' Number antibodies gamma 1 '''
sb.set(context='paper', style='ticks')
g = sb.FacetGrid(data=data, col='cov', hue='num_abs', palette='RdPu', height=height, aspect=aspect, legend_out=True,col_wrap=4)
g.map(plt.plot,labels[0], labels[1], marker='s', ms=ms,ls='', alpha=0.7)
plt.ylim([1,ymax])
plt.semilogy()
plt.semilogx()
plt.tight_layout() 
plt.savefig('../../output/figs/cocktail_num_abs.png', dpi=300, transparent=True)
plt.show()

''' Infective virus coverage gamma 2 '''
sb.set(context='paper', style='ticks')
g = sb.FacetGrid(data=data, col='cov', hue='inf_cov', palette='OrRd', height=height, aspect=aspect, legend_out=True,col_wrap=4)
g.map(plt.plot,labels[0],labels[1], marker='s', ms=ms,ls='', alpha=0.7)
plt.ylim([1,ymax])
plt.semilogy()
plt.semilogx()
plt.tight_layout() 
plt.savefig('../../output/figs/cocktail_inf_cov.png', dpi=300, transparent=True)
plt.show() 

''' Virus coverage gamma 2 '''
sb.set(context='paper', style='ticks')
g = sb.FacetGrid(data=data, col='cov', hue='opt_cov', palette='YlGn', height=height, aspect=aspect, legend_out=True,col_wrap=4)
g.map(plt.plot,labels[0],labels[1], marker='s', ms=ms,ls='', alpha=0.7)
plt.ylim([1,ymax])
plt.semilogy()
plt.semilogx()
plt.tight_layout() 
plt.savefig('../../output/figs/cocktail_opt_cov.png', dpi=300, transparent=True)
plt.show() 

hue_order = data.cocktail.unique()

'''Cocktail names ''' 
sb.set(context='paper', style='ticks')
g = sb.FacetGrid(data=data, col='cov', hue='cocktail', hue_order = hue_order , palette='Set3', height=height, aspect=aspect, legend_out=True,col_wrap=4)
g.map(plt.plot,labels[0], labels[1], marker='s', ms=ms,ls='', alpha=0.7)
plt.ylim([1,ymax])
plt.semilogy()
plt.semilogx()
#g.add_legend()
plt.tight_layout() 
plt.savefig('../../output/figs/cocktail_names.png', dpi=300, transparent=True)
plt.show()


print(data.head())

'''Is optimal cocktail ''' 
sb.set(context='paper', style='ticks')
g = sb.FacetGrid(data=data, col='cov', hue='isopt' , palette='Blues', height=height, aspect=aspect, legend_out=True,col_wrap=4)
g.map(plt.plot,labels[0], labels[1], marker='s', ms=ms,ls='', alpha=0.7)
plt.ylim([1,ymax])
plt.semilogy()
plt.semilogx()
#g.add_legend()
plt.tight_layout() 
plt.savefig('../../output/figs/cocktail_isopt.png', dpi=300, transparent=True)
plt.show()



''' Get antibody mix properties ''' 
cocktail_file= '../../output/cocktail/cocktails_'
cocktail_files = [cocktail_file + str(cov)+'.csv' for cov in coverage]


#cocktail.graph_robustness()
#exit()
#cocktail.get_antibody_mix_properties(cocktail_files)

#cocktail.graph_cocktail_mixes()

#cocktail.graph_cocktail_distance()
#exit()
#cocktail.graph_UMAP_cocktail()

#exit()
''' Graph individual simulations ''' 
i = 10
'''
for i in range(len(coverage)):
    cocktail_file= '../../output/cocktail/cocktails_' + str(coverage[i]) +'.csv'
    cocktail.graph_simulation(cocktail_file)
'''

#cocktail_file= '../../output/cocktail/cocktails_'

''' Graph gamma 1 vs number of antibodies '''
#cocktail_files = [ cocktail_file + str(cov) +'.csv' for cov in coverage]
#cocktail.graph_optimization_results(cocktail_files, gamma1_range=None, gamma2_range=None)

#exit()
''' Graph heatmap of ab cocktail mixes '''

#cocktail.graph_cocktail_mixes()
#exit()


''' Get maximum virus coverage, get maximum neutralization coverage  for each simulation '''
#cocktail.compare_cocktail_sensitivities(cocktail_file, c1, c2)


#dfm = dfc.melt(id_vars='mut', value_vars=dfc.columns[1:])
#dfm['location'] = dfm.mut.str[1:-1]
#dfm = dfm.rename(columns={'variable':'antibody', 'value':'ddg'})


#classes: C105_TH28I_YH58F'class I , 'C135'class III, 'REGN10987' class III, 'COVA2-39' class II
#colors = ['#80b1d3','#fc8d62', '#66c2a5','#ffd92f' ]
#hue_order  = ['C105_TH28I_YH58F', 'C135', 'COVA2-39','REGN10987']

#plt.figure(figsize=(2.5,2.5))
#sb.set(context='paper', style='white')
#ax = sb.scatterplot(data=merged, y='ddgn', x = 'ACE2', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.5)
#plt.axhline(y=0, color='grey',xmin=-1, xmax=1, lw=0.5)
#plt.axvline(x=0, color='grey',ymin=-1, ymax=1, lw=0.5)
#ax.legend().remove()
#plt.xticks(rotation=0)
#plt.tight_layout() 
#plt.savefig('../../output/manuscript/figs/combination_sensitivity.png', dpi=300)
#plt.show()    

''' Graph ACE2 ddg greater than zero ''' 


''' Grouped by count  '''
#grouped = dfmp.groupby('antibody').count().reset_index()
#total_muts = len(dfm.mut.unique())
#grouped['pc']  = (total_muts -grouped.ddg)/total_muts*100.

''' Grouped by mean antibody at location ''' 
#grouped_mean = dfmp.groupby(['location', 'antibody']).mean().reset_index()
#print(grouped_mean.head())

grouped_mean.to_csv('../../output/tmp.csv')

plt.figure(figsize=(5.5,2.5))
sb.set(context='paper', style='white')
#ax = sb.stripplot(data=dfmp, x='location', y = 'ddg', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.3, jitter=True)
ax = sb.stripplot(data=dfmp, x='mut', y = 'ddg', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.3, jitter=True)
ax.legend().remove()
plt.xticks(rotation=90)
plt.tight_layout() 
plt.savefig('../../output/manuscript/figs/combination_sensitivity.png', dpi=300)
plt.show()    
exit()

sb.set(context='paper', style='ticks')
sb.violinplot(data=dfmp, x='location',y = 'ddg', color ='white')
ax=sb.stripplot(data=dfmp, x='location',y = 'ddg',alpha = 0.5, hue='antibody', palette=colors)
plt.xticks(rotation=90)
ax.legend().remove()
plt.tight_layout()
plt.show()

plt.figure(figsize=(1.5,1.25))
sb.set(context='paper', style='ticks')
sb.barplot(data=grouped, x='antibody',y = 'pc', palette=colors)
#plt.xticks([])
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('../../output/manuscript/figs/vcov.png', dpi=300)
plt.show()

 
'''
with open('../../output/manuscript/figs/opt_results.txt', 'w') as writer:
    for choice in choices:
        writer.write(choice+'\n')
'''


exit()


''' Analysis post running of algorithm ''' 



exit()


