import numpy as np
import cvxpy as cp
import pandas as pd
import seaborn as sb 
import matplotlib.pyplot as plt

def compute_antibodies(fitness_matrix, k, lmbda1, lmbda2, virus_matrix):
    m,n = fitness_matrix.shape
    num_unmanaged = int(m*k)
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix).clip(min=0)
    virus_incidence_matrix = np.sign(virus_matrix).clip(min=0)
    constraints = [cp.sum(c) >= 1, cp.sum_smallest((incidence_matrix@c), num_unmanaged+1) >= 1] # c[0] == 1 ==> add to list if want first antibody chosen
    objective = cp.Minimize(lmbda1*cp.norm1(c)+cp.matmul(cp.sum(fitness_matrix, axis=0), c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value


def import_fitness_landscapes():
    filepath = '../data/'
    combinepath = '../output/combine/'
    
    filename = 'rbd_ab_fitness'
    virus_filename = 'rbd_ace2_fitness'


    df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0, index_col=0)
    virus_df = pd.read_csv(filepath+virus_filename+'.csv', sep=',', header=0, index_col=0)
    virus_df = virus_df[~df.isnull().any(axis=1).values]
    df = df.dropna()
    virus_matrix = virus_df.values
    fitness_matrix = df.values

    return df, virus_df, fitness_matrix, virus_matrix 


''' Run all the simulations '''
def run_simulations():


    df,vdf, fitness_matrix, virus_matrix = import_fitness_landscapes()

    allresults = pd.DataFrame()
    allresults['antibody'] = df.columns

    lmbda1 = 580
    lmbda2 = 1

    ks = [0.06,0.08,0.09,0.1,0.15,0.20,0.30,0.4,0.5, 0.4, 0.3]
    for k in ks:
        results = compute_antibodies(fitness_matrix, k, lmbda1, lmbda2, virus_matrix)
        choices = df.columns[results.astype(bool)]
        allresults[str(1-k)]=results
        
        allresults.to_csv('../output/combine/allresults.csv', index=False)

''' Graph simulations ''' 
def graph_simulations():

    allres = pd.read_csv('../output/combine/allresults.csv')
    
    print(allres.head())
    allres = allres.replace({'C105_TH28I_YH58F':'C105$^{T28IY58F}$'})
    melted = allres.melt(id_vars='antibody', value_vars= allres.columns[1:])
    allres= allres.set_index('antibody')
    allres = allres.drop('C002-S2', axis=0)
    allres = allres.drop('0.5', axis=1)

    print(allres.columns)
    cols = allres.columns.astype(float)*100
    
    cols = cols.astype(int).astype(str) +'%'
    
    dictcols = dict(zip(allres.columns, cols))
    allres = allres.rename(columns=dictcols)
    print (cols) 
    allres = allres.transpose()

    figsize = (2.5,5.5)
    cmap = ['#fde0dd','#c51b8a']
    plt.figure(figsize= figsize)
    sb.set(context='paper', style='white')
    sb.heatmap(data=allres ,cmap=cmap, square=True, linecolor='white', linewidths=0.5, cbar=None, alpha=0.7)
    plt.xticks(rotation = 90)
    plt.tight_layout() 
    plt.savefig('../output/figs/combinations.png', dpi=300)
    plt.show()    


def analyze_cocktail(df,vdf, cocktail, cocktailname, graph=True):
    
    ''' Bar graph with the antibodies that are chosen virus coverage '''
    dfc = df[cocktail]
    dfc = dfc.reset_index()
    
    ''' Find min in rows '''
    dfc['min_ddg'] = dfc.iloc[:,1:].min(axis=1)
    dfc['ispos'] = dfc.min_ddg>0 
    
    cov = float(cocktailname)
    allmuts = len(dfc.mut.unique())
    uncovmuts = len(dfc.loc[dfc.ispos].mut.unique())
    covmuts = allmuts - uncovmuts 
    descov = allmuts*(cov)

    dfc['location'] = dfc.mut.str[1:-1]
    
    ''' All positive ispos ddgs are uncovered''' 
    dfu = dfc.loc[dfc.ispos]
    dfu['ddgn'] = dfu.min_ddg/dfu.min_ddg.max()    
    
    merged = vdf.merge(dfu, on='mut', how='inner')
    merged['minACE2']= -merged.ACE2
    merged['cocktail'] = cocktailname 

    ''' Do scatter for each antibody '''
    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='white')
    ax = sb.scatterplot(data=merged, y='ddgn', x = 'minACE2', alpha=0.8,hue ='location', palette='Set3')
    plt.axhline(y=0, color='grey',xmin=-1, xmax=1, lw=0.5)
    plt.axvline(x=0, color='grey',ymin=-1, ymax=1, lw=0.5)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.ylabel('ddg(Ab/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_sensitivity' + cocktailname + '.png', dpi=300)

    ''' Do histogram/KDE plot '''
    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='white')
    ax = sb.histplot(data=merged,  x = 'minACE2', alpha=0.5)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_kde' + cocktailname + '.png', dpi=300)

    ''' Violin plot of all sensitive locations '''
    plt.figure(figsize=(3.5,1.5))
    sb.set(context='paper', style='ticks')
    sb.violinplot(data=merged, y='minACE2', x = 'location', color='white',scale='width', inner='point')
    ax = sb.stripplot(data=merged, y='minACE2', x = 'location', alpha=0.5, color='Red', jitter=True)
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_infectionsensitivity_' + cocktailname +'.png', dpi=300)

    if graph: 
        plt.show()
        
    return merged 


def compare_cocktail_sensitivities():

    allres = pd.read_csv('../output/combine/allresults.csv')

    ''' Pick ones to compare ''' 
    c1name = '0.85'
    c2name = '0.6'
    allres = allres[['antibody',c1name, c2name]]

    c1 = allres.loc[allres[c1name] == 1].antibody.values 
    c2 = allres.loc[allres[c2name] == 1].antibody.values 

    print(c1)
    print(c2)

    ''' graph virus sensitivities '''
    df,vdf, fitness_matrix, virus_matrix = import_fitness_landscapes()
    c1data = analyze_cocktail(df, vdf,c1, c1name, graph=False)
    c2data = analyze_cocktail(df, vdf,c2, c2name, graph=False)


    merged = pd.concat([c1data, c2data], axis=0)
    print(c1data.head())
    print(c2data.head())

    colors = ['#de77ae', '#7fbc41']
    hue_colors = ['0.85','0.6']
    ''' Do histogram/KDE plot '''
    plt.figure(figsize=(1.75,1.75))
    sb.set(context='paper', style='ticks')
    ax = sb.histplot(data=merged,  x = 'minACE2', alpha=1, hue='cocktail', palette= colors, hue_order=hue_colors)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_kde.png', dpi=300)
    
    print(merged.head())
    sorted= merged.sort_values('location')

    print(sorted.head())

    plt.figure(figsize=(3.5,1.5))
    sb.set(context='paper', style='white')
    ax =sb.violinplot(data=sorted, y='minACE2', x = 'location', color='white',scale='width', inner='point', hue='cocktail', hue_order=hue_colors,palette=colors, legend=False, lw=0.5,alpha=0.7)
    #ax = sb.scatterplot(data=sorted, y='minACE2', x = 'location', alpha=0.5, palette=colors, hue='cocktail')
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_sensitivity.png', dpi=300)


    ''' Group the cocktails '''
    cmean = sorted.groupby(['location','cocktail']).mean().reset_index()[['cocktail','location','minACE2','ddgn']]
    cmean = cmean.loc[cmean.minACE2 <=0]
    sorted = sorted.loc[sorted.minACE2 <= 0]



    print(cmean.location.unique())
    plt.figure(figsize=(1.5,2.5))
    sb.set(context='paper', style='white')
    #ax =sb.violinplot(data=sorted, y='minACE2', x = 'location', color='white',scale='width', inner='point', hue='cocktail', palette='Set2', legend=False, lw=0.5,alpha=0.7)
    #ax = sb.scatterplot(data=cmean, y='minACE2', x = 'location', alpha=0.8, palette=colors, hue='cocktail',hue_order=hue_colors, size='ddgn')
    ax = sb.scatterplot(data=sorted, x='minACE2', y = 'location', alpha=0.8, palette=colors, hue='cocktail',hue_order=hue_colors, s=30, marker='s', size='min_ddg')
    #ax = sb.stripplot(data=sorted, y='minACE2', x = 'location', alpha=0.7, palette=colors, hue='cocktail',hue_order=hue_colors, dodge=True)
    #ax = sb.displot(data=sorted, y='minACE2', x = 'location', alpha=0.8, palette=colors, hue='cocktail',hue_order=hue_colors)
    #ax = sb.barplot(data=sorted, y='minACE2', x = 'location', alpha=0.5, palette=colors, hue='cocktail',hue_order=hue_colors)
    #ax = sb.scatterplot(data=cmean, y='minACE2', x = 'ddgn', alpha=0.5, palette=colors, hue='cocktail',hue_order=hue_colors)

    #ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    #plt.xticks(rotation=90)
    #plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_sensitivity.png', dpi=300)
    plt.show()

    print(cmean.head())
    exit()


    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='ticks')
    #ax =sb.violinplot(data=sorted, y='minACE2', x = 'location', color='white',scale='width', inner='point', hue='cocktail', palette='Set2', legend=False, lw=0.5,alpha=0.7)
    ax = sb.kdeplot(data=sorted, x='minACE2', alpha=0.5, palette=colors, hue='cocktail')
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../output/figs/combination_scatter.png', dpi=300)



#graph_simulations()

compare_cocktail_sensitivities()
exit()


dfm = dfc.melt(id_vars='mut', value_vars=dfc.columns[1:])
dfm['location'] = dfm.mut.str[1:-1]
dfm = dfm.rename(columns={'variable':'antibody', 'value':'ddg'})


#classes: C105_TH28I_YH58F'class I , 'C135'class III, 'REGN10987' class III, 'COVA2-39' class II
colors = ['#80b1d3','#fc8d62', '#66c2a5','#ffd92f' ]
hue_order  = ['C105_TH28I_YH58F', 'C135', 'COVA2-39','REGN10987']

plt.figure(figsize=(2.5,2.5))
sb.set(context='paper', style='white')
ax = sb.scatterplot(data=merged, y='ddgn', x = 'ACE2', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.5)
plt.axhline(y=0, color='grey',xmin=-1, xmax=1, lw=0.5)
plt.axvline(x=0, color='grey',ymin=-1, ymax=1, lw=0.5)
ax.legend().remove()
plt.xticks(rotation=0)
plt.tight_layout() 
plt.savefig('../output/figs/combination_sensitivity.png', dpi=300)
plt.show()    

''' Graph ACE2 ddg greater than zero ''' 




''' Grouped by count  '''
grouped = dfmp.groupby('antibody').count().reset_index()
total_muts = len(dfm.mut.unique())
grouped['pc']  = (total_muts -grouped.ddg)/total_muts*100.

''' Grouped by mean antibody at location ''' 
grouped_mean = dfmp.groupby(['location', 'antibody']).mean().reset_index()
print(grouped_mean.head())

grouped_mean.to_csv('../output/tmp.csv')

plt.figure(figsize=(5.5,2.5))
sb.set(context='paper', style='white')
#ax = sb.stripplot(data=dfmp, x='location', y = 'ddg', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.3, jitter=True)
ax = sb.stripplot(data=dfmp, x='mut', y = 'ddg', hue ='antibody', hue_order=hue_order, palette=colors, alpha=0.3, jitter=True)
ax.legend().remove()
plt.xticks(rotation=90)
plt.tight_layout() 
plt.savefig('../output/figs/combination_sensitivity.png', dpi=300)
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
plt.savefig('../output/figs/vcov.png', dpi=300)
plt.show()

 

with open(combinepath+'opt_results.txt', 'w') as writer:
    for choice in choices:
        writer.write(choice+'\n')
