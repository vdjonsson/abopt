import random 
import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
import plotutils as pu 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import gcf
from sklearn.manifold import MDS, TSNE
import numpy as np 
import scipy.stats as stats 
from sklearn import preprocessing
import logomaker
import manuscript_analysis as ma
from scipy import stats

''' This file will run all the figure plotting for the paper  ''' 

''' This is a biorender image ''' 
def figure1a(): 
    print('This is a biorender image') 


''' Natalie: add pointer to run your code here ''' 
def figure1b():
    ''' Natalie: add pointer to run your code here ''' 
    print ('ND add code') 

''' logoplots '''
def figure1c():

    crp_df = -logomaker.get_example_matrix('crp_energy_matrix',
                                        print_description=False)
    print(crp_df)
    # create Logo object
    crp_logo = logomaker.Logo(crp_df,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold')

    # style using Logo methods
    crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    
    # style using Axes methods
    crp_logo.ax.set_ylabel("$-\Delta \Delta G$ (kcal/mol)", labelpad=-1)
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.show()


def figure_2c():

    alg = 'MDS'
    #alg = 'UMAP'
    print('figure 2c') 

    ''' Read ab configuration file '''
    abdf = pd.read_table('../../manuscript/antibody_list.txt', sep=',') 
    antibody_virus_landscape = pd.read_csv('../../output/merge/rbd_ab_fitness.csv')    

    abt = antibody_virus_landscape.transpose() 
    abs = antibody_virus_landscape.columns[1:].values
    print(abs)

    abdf = abdf.loc[abdf.antibody.isin(abs)]
    print(abdf.head())
    ab_class = abdf['abclass'].values 
    ab_names = abdf['antibody'].values 
    ab_vhgene = abdf['vhgene'].values 

    ab_class_dict = dict(zip(ab_names, ab_class))
    ab_gene_dict = dict(zip(ab_names, ab_vhgene))

    classes = [ab_class_dict[x] for x in abs]
    genes = [ab_gene_dict[x] for x in abs]

    ''' First normalize the data ''' 

    ''' Graph embedding antibody dimensionality reduction '''
    embedding = ma.learn_antibody_distance(antibody_virus_landscape, alg = alg)

    embedding['ab'] = abs
    embedding['ab_class'] = classes
    embedding['ab_gene'] = genes
    print (embedding)
    
    class_color = sb.color_palette('colorblind', 4)
    colors = class_color 
    hue = 'ab_class'

    ''' Plot ''' 
    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(4.2,4))
    #f,ax = plt.subplots()
    sb.set(context='paper')
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue, palette=colors, s=80, legend=False)
    sb.despine()
    plt.gca().set_aspect('equal', 'datalim')
    for j, lab in enumerate(embedding.ab):
        ax.annotate(lab, (embedding.dim1[j]-1.5, embedding.dim2[j]+ 1.5 ))                
    plt.title(alg)
    plt.tight_layout()
    plt.savefig('../../output/figs/' + alg + '_virus_scan_WT_ab_class.png', dpi = 300)
    plt.show()

    gene_color = sb.color_palette('Set2', 6)
    colors = gene_color 
    hue = 'ab_gene'
    hue_order = ['VH3-53', 'VH3-30', 'VH1-46', 'VH1-2','VH5-51', 'Unknown'] 

    ''' Plot ''' 
    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(4.2,4))
    #f,ax = plt.subplots(figsize=(6,6))
    #f,ax = plt.subplots()
    sb.set(context='paper')
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue,hue_order=hue_order, palette=colors, s=80, legend='full',facecolor='white')
    plt.legend(facecolor='white')
    sb.despine()
    plt.gca().set_aspect('equal', 'datalim')
    for j, lab in enumerate(embedding.ab):
        ax.annotate(lab, (embedding.dim1[j]+1.5, embedding.dim2[j]+ 1.5 ))                
    plt.title(alg)
    plt.tight_layout()
    plt.savefig('../../output/figs/' + alg + '_virus_scan_WT_ab_gene.png', dpi = 300)
    plt.show()


def figure_4b():

    rbd = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_with_predictors.csv')

    ''' Read RBD fitness landscape ''' 
    rbd = rbd[['site_SARS2', 'bind_avg']]
    rbd['-bind_avg'] = -1 * rbd.bind_avg.values  
    vmin  = rbd['-bind_avg'].min()
    vmax  = rbd['-bind_avg'].max()

    ''' Read random forest features '''
    rf  = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_filtered.csv')
    featurelocs = rf.feature.unique()

    ''' Calculate mean ''' 
    grouped = rbd.groupby('site_SARS2').mean().reset_index()
    grouped_small = grouped.loc[grouped.site_SARS2.isin(featurelocs)]

    ''' Get epitope locations '''
    eps = pd.read_csv('../../data/location/epitopes.csv').epitope_location.values
    eps_rbd = []
    for loc in grouped_small.site_SARS2: 
        isep = False
        if loc in eps: 
            isep = True
        eps_rbd.append(isep)

    eps_labels = pd.DataFrame(eps_rbd,index=grouped_small.site_SARS2).transpose()

    print(eps_labels)

    cmap = [ '#feb24c','#f03b20']
    plt.figure(figsize=(10,1))
    sb.set(context='paper', style='ticks') 
    sb.heatmap(data=eps_labels, cmap=cmap)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('../../output/figs/' + 'figure_4b_ace2'+ '.png', dpi=300, transparency=True)
    plt.show()
    
    
    ''' Plot '''
    hm = grouped.set_index('site_SARS2')
    hm = hm[['-bind_avg']].transpose()

    hm_small = grouped_small.set_index('site_SARS2')
    hm_small = hm_small[['-bind_avg']].transpose()

    plt.figure(figsize=(10,1))
    sb.set(context='paper', style='ticks') 
    sb.heatmap(data=hm, cmap='Reds',vmin=vmin , vmax=vmax)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10,4))
    sb.set(context='paper', style='ticks') 
    sb.heatmap(data=hm_small, cmap='RdBu_r',vmin=vmin , vmax=vmax)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('../../output/figs/' + 'figure_4b'+ '.png', dpi=300, transparency=True)
    plt.show()

    ''' Plot '''
    sb.set(context='paper', style='ticks') 
    sb.violinplot(data=rbd, y='-bind_avg', x='site_SARS2', color='Red')
    plt.xticks(rotation=90)
    plt.show()

    print(merged.head())


def calculate_C105_optimized_landscape(data):
    ''' C105 optimized epitope  '''
    epi = '6xcm_Repair_TH28I_YH58F_epitopes.csv'
    epdf  = pd.read_csv('../../output/epitope/'+ epi)

    epitope = epdf['number_virus'].unique()
    muts = data.loc[data['location'].isin(epitope)].mut.values

    data['C105_TH28I_YH58F_new'] = data.C105.values

    for mut in muts: 

        new_val = data.loc[data.mut==mut].C105_TH28I_YH58F.values[0]
        old_val = data.loc[data.mut==mut].C105.values[0]
        data.loc[data['mut'] == mut, 'C105_TH28I_YH58F_new'] = new_val - old_val

    data = data.drop('C105_TH28I_YH58F', axis=1)
    data = data.rename(columns={'C105_TH28I_YH58F_new':'C105_TH28I_YH58F'})

    return data

 
def clean_filter_landscape(data): 

    locations = data.mut.str[1:-1].astype(int)
    data['location'] = data.mut.str[1:-1].astype(int)
    data = data.loc[data['location'] <520]
    ab_order = ['mut','location','C105','C105_TH28I_YH58F','B38','CB6','CV30','CC121','C002-S1', 'C002-S2','C119','C121-S1','C144','COVA2-39','C135','C110' ,'REGN10987','REGN10933', 'ACE2']
    data = data[ab_order]
    return data 

def normalize_landscape(data):
    
    ''' Normalize antibodies '''
    ''' [-1 1] every column '''

    data = data.fillna(0)
    values = data.iloc[:,2:].values
    print(data.columns)
    
    ''' Normalize l2 this should be changed '''
    normalizer = preprocessing.Normalizer()
    data_norm = pd.DataFrame(normalizer.transform(values), columns=data.columns[2:], index=data['location'])
    data_norm['mut'] = data.mut.values

    return data_norm
    
def get_antibody_fitness_landscape(): 

    merge_dir = '../../output/merge/'
    data = pd.read_csv(merge_dir + 'rbd_ace2_ab_fitness.csv')
    ab_list = pd.read_table('../../manuscript/antibody_list.txt', sep=',')

    data_clean = clean_filter_landscape(data) 
    data_opt = calculate_C105_optimized_landscape(data_clean)
    data_norm = normalize_landscape(data_opt)
    
    data_norm = data_norm.reset_index()
    return data_norm


def figure_4c(): 

    df = get_antibody_fitness_landscape()
    df = df.drop('ACE2', axis=1)

    ''' Calculate mean ''' 
    grouped = df.groupby('location').mean().reset_index()
    grouped = grouped.set_index('location').transpose()

    mask = np.zeros_like(grouped)
    mask = grouped == 0  # mean equal to zero signifies that the structure is missing

    ab_list = pd.read_table('../../manuscript/antibody_list.txt', sep=',')

    ''' Get all the class and VH gene mappings ''' 
    ab_vh = dict(zip(ab_list.antibody.values, ab_list.vhgene.values))
    ab_class = dict(zip(ab_list.antibody.values, ab_list['abclass'].values))
    ab_vh_gene = [ab_vh[ab] for ab in grouped.index]
    ab_class = [ab_class[ab] for ab in grouped.index]

    grouped['vh']= ab_vh_gene
    grouped['abclass']= ab_class

    vhgene = grouped.pop('vh')
    abclass = grouped.pop('abclass')

    pal = sb.color_palette('Set2', 6)
    pal_class = sb.color_palette('colorblind', 4)
    lut = dict(zip(vhgene.unique(), pal))
    lut_class = dict(zip(abclass.unique(), pal_class))
    vh_colors = vhgene.map(lut)
    class_colors = abclass.map(lut_class)

    row_colors = pd.DataFrame()
    row_colors['Class'] = class_colors 
    row_colors['VH'] = vh_colors 

    ''' Get RBD overlap with ACE2, get overlap with features'''
    fs = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_with_predictors.csv')
    fs2 = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_coefficients.csv')
    fs['wt'] = fs.mutation.str[0] == fs.mutation.str[-1] 
    fss = fs.loc[fs.wt == True]

    print(fss.head())
    fss['wt_num'] = fss.mutation.str[0:-1]
    print(fss.head())

    ''' Epitope location ''' 
    '''
    dir = '../../data/pdb/RBD_ACE2/'
    fname = '6m0j.pdb'
    distance = 6.0

    labeled_chains = structure.label_chains(fname[0:-4]) 
    print(labeled_chains)
    epitopes = structure.find_epitopes(dir, fname, labeled_chains, distance)
    '''

    eps = pd.read_csv('../../data/location/epitopes.csv').epitope_location.values
    eps_rbd = []
    for loc in grouped.columns: 
        isep = False 
        if loc in eps: 
            isep = True 
        eps_rbd.append(isep)

    print (eps_rbd) 
    eps_labels = pd.Series(eps_rbd,index=grouped.columns)

    eps_pal = sb.color_palette('YlOrRd', eps_labels.unique().size)
    eps_lut = dict(zip(eps_labels.unique(), eps_pal))
    eps_colors = eps_labels.map(eps_lut)
    
    rbd_ace = fs2.merge(fss, how='inner', right_on='wt_num', left_on='rbd_mutation')
    rbd_ace['location'] = rbd_ace.rbd_mutation.str[1:].astype(int)
    rbd_ace = rbd_ace[['location', 'importances']]
    rbd_ace = rbd_ace.set_index('location')
    rbd_ace['RF'] = rbd_ace.importances > 1e-2
    rbd_ace = pd.Series(rbd_ace.RF.astype(str), index= rbd_ace.index)

    rbd_labels = rbd_ace 
    rbd_pal = sb.color_palette('YlOrRd', rbd_labels.unique().size)
    rbd_lut = dict(zip(rbd_labels.unique(), rbd_pal))
    rbd_colors = rbd_labels.map(rbd_lut)

    col_colors = pd.DataFrame()
    col_colors = eps_colors
    col_colors = pd.DataFrame()
    col_colors['ACE2 epitope'] = eps_colors 

    grouped.to_csv('../../output/tmp.csv') 

    sb.set(context='paper', font_scale=0.65)
    figsize= (8,2.75)
    g = sb.clustermap(data=grouped, cmap='RdBu_r', vmin=-1, vmax=1, row_cluster = True, col_cluster=False,figsize=figsize, row_colors=row_colors, col_colors=col_colors, mask=mask)
    
    plt.title('ddg(ab/RBD)<0 normalized')

    for label in abclass.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut_class[label], label=label, linewidth=0)

    l1 = g.ax_col_dendrogram.legend(title='', loc="center", ncol=12, bbox_to_anchor=(0.5, 0.9), bbox_transform=gcf().transFigure, facecolor ='white')

    for label in vhgene.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)

    l2 = g.ax_col_dendrogram.legend(title='', loc="center", ncol=6, bbox_to_anchor=(0.8, 0.9), bbox_transform=gcf().transFigure, facecolor ='white')

    '''
    for label in rbd_labels.unique():
        g.ax_col_dendrogram.bar(2, 0, color=rbd_lut[label], label=label, linewidth=0)

    l2 = g.ax_col_dendrogram.legend(title='class/VH', loc="center", ncol=6, bbox_to_anchor=(0.5, 0.9), bbox_transform=gcf().transFigure)
    '''
    for label in eps_labels.unique():
        g.ax_col_dendrogram.bar(0, 0, color=eps_lut[label], label=label, linewidth=0)

    l3 = g.ax_col_dendrogram.legend(title='', loc="center", ncol=6, bbox_to_anchor=(0.5, 0.9), bbox_transform=gcf().transFigure, facecolor='white')
    
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig('../../output/figs/' + 'figure_4d'+ '.png', dpi=300, transparency=True)
    plt.show()

    exit()
    ''' Calculate number of positive ddgs '''
    dfs = dfs > 0
    #dfs['mut'] = df.mut
    dfs['location'] = df.index

    print(dfs.head())
    grouped = dfs.groupby(by='location').sum().reset_index()
    grouped = grouped.set_index('location').transpose()

    sb.set(context='paper')
    plt.figure(figsize=(12,2.5))
    g = sb.heatmap(data=grouped, mask=mask)#, cmap ='RdBu_r')
    g.set_facecolor('#756bb1')
    plt.tight_layout()
    plt.title('ddg(ab/RBD)<0 normalized')
    plt.tight_layout()
    plt.savefig('../../output/figs/'+ 'hm_count'+ '.png', dpi=300)
    #plt.show()

def get_antibody_properties():

    abdf = pd.read_table('../../manuscript/antibody_list.txt', sep=',') 
    ab_class = abdf['abclass'].values 
    ab_names = abdf['antibody'].values 
    ab_genes = abdf['vhgene'].values 

    ab_class_dict = dict(zip(ab_names, ab_class))
    ab_gene_dict = dict(zip(ab_names, ab_genes))

    return ab_class_dict, ab_gene_dict 

def melt_landscape(data):

    melted = data.melt(value_vars=data.columns[1:-1], id_vars=['location', 'mut'])
    melted = melted.rename(columns={'variable':'antibody', 'value':'binding'})
    return melted

def figure_4e(): 

    data = get_antibody_fitness_landscape()
    melted = melt_landscape(data)

    classes, genes = get_antibody_properties()
 
    for ab in list(genes): 
        melted.loc[melted['antibody'] == ab, 'vhgene'] = genes[ab] 
        melted.loc[melted['antibody'] == ab, 'class'] =  classes[ab] 


    print(melted.head())

    
    ''' Reduce set '''
    filtered = melted.loc[melted.antibody.isin(['C105_TH28I_YH58F','C105'])]


    ''' Perform t test and compare accross locations designed ab vs wt  '''
    ttest_vector = []
    stat_test = pd.DataFrame()

    for loc in filtered['location'].unique():

        filtered_loc = filtered.loc[filtered['location'] == loc]
        opt = filtered_loc.loc[filtered_loc.antibody == 'C105_TH28I_YH58F'].binding
        old = filtered_loc.loc[filtered_loc.antibody == 'C105'].binding

        res = stats.ttest_ind(opt, old)  
        stat_test = stat_test.append(pd.Series([loc, res.pvalue]), ignore_index=True)



    ''' Statistically significant locations to compare '''
    stat_test = stat_test.rename(columns={0:'location', 1:'pval'})
    stat_test = stat_test.loc[stat_test.pval < 5e-2]
    difflocs = stat_test.location.astype(int).values

    print (difflocs)

    ''' Plot locations of epitope only to compare '''
    epi = '6xcm_Repair_TH28I_YH58F_epitopes.csv'
    epdf  = pd.read_csv('../../output/epitope/'+ epi)
    epitopes = epdf['number_virus'].unique()

    ''' Significant locations '''
    mystr = 'locs'
    dfd = filtered.loc[filtered['location'].isin(epitopes)]
    dfd = filtered.loc[filtered['location'].isin(difflocs)]

    sig = [x in difflocs for x in dfd['location'].values]

    dfd['significance'] = sig
    sigloc = dfd['location'].unique() 

    colors = ['#df65b0', '#74a9cf']


    ''' Plot locations against ace2 binding '''
    fit = pd.read_csv('../../output/merge/rbd_ace2_ab_fitness.csv')[['ACE2','mut']]
    
    merged = data.merge(fit, how='inner')

    grouped = merged.groupby('location').mean().reset_index()
    grouped = merged
    print(grouped) 
    
    col1 ='C105_TH28I_YH58F'
    col2 ='C105'

    ''' Just use epitopes location ''' 
    grouped = grouped.loc[grouped['location'].isin(epitopes)]

    
    ''' Plot scatter'''
    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='ticks', font_scale=1.3) 
    ax = sb.scatterplot(data=grouped, x=col2, y=col1,color='Blue', hue='location')
    #ax = sb.scatterplot(data=grouped, x='ACE2', y=col2,color ='Red' , hue ='location', marker ='D' ,)
    #ax = sb.scatterplot(data=merged, x='ACE2', y='C105_TH28I_YH58', palette=colors)
    ax.set_alpha(0.3)
    plt.axhline(0, lw=0.5, color = 'gray')
    plt.axvline(0, lw=0.5, color = 'gray')
    #plt.xlim([-1,1])
    #plt.ylim([-1,1])
    ax.legend().set_visible(True)
    plt.tight_layout()   
    plt.savefig('../../output/figs/' + 'design_fitland_virus_scan_WT_ab_bind' + mystr +'.png', dpi = 300, transparent=True)
    plt.show()

    
    ''' Plot locations '''
    plt.figure(figsize=(6.5,2.5))
    sb.set(context='paper', style='ticks', font_scale=1.3) 
    ax = sb.violinplot(data=dfd, x='location', y='binding',alpha=0.1,scale='width', hue='antibody',palette=colors, legend=False, split=True, inner='stick')
    ax.set_alpha(0.3)
    plt.axhline(0, lw=0.5, color = 'gray')
    plt.xticks(rotation=90)
    ax.legend().set_visible(True)
    plt.tight_layout()   
    plt.savefig('../../output/figs/' + 'design_fitland_virus_scan_WT_ab_bind' + mystr +'.png', dpi = 300, transparent=True)
    plt.show()

    property = 'vhgene'
    palette = 'Set2'

    #property = 'class'
    #palette = 'colorblind'

    #property = 'antibody'
    #palette = 'colorblind'

    print(melted.head())
    melted = melted.loc[melted.binding> 0]

    g = sb.FacetGrid(melted, row=property,aspect=8.5, height=0.8, palette=palette, hue=property)        
    #g.map_dataframe(sb.kdeplot, data= df, x ='location', clip_on=False,fill=False, alpha=0.8, lw=1.5)
    g.map_dataframe(sb.histplot, data= melted, x ='location',stat='density',fill=True, alpha=0.6,bins=100)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, color=color, ha="left", va="center", transform=ax.transAxes)

    g.map(label,'antibody')
    
    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-10)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    plt.tight_layout()
    plt.savefig('../../output/figs/' + 'fitland_virus_scan_WT_ab_' + property +'.png', dpi = 300)
    plt.show()
    exit()

    sb.set(context='paper', style='ticks')
    sb.kdeplot(data=data, x='location',hue = property, fill=True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    sb.set(context='paper', style='ticks')
    g = sb.FacetGrid(data=data, col="class", col_wrap=2, legend_out=True) # , height=1.5, aspect=.5)
    g.map_dataframe(sb.kdeplot, x='location',hue = property,shade=True)
    g.add_legend()
    plt.tight_layout()
    plt.show()

    sb.set(context='paper', style='ticks')
    g = sb.FacetGrid(data=data, col="vhgene", col_wrap=2, legend_out=False) # , height=1.5, aspect=.5)
    g.map_dataframe(sb.kdeplot, x='location',hue = 'class',shade=True).add_legend()
    plt.tight_layout()
    plt.show()

    sb.set(context='paper', style='ticks')
    g = sb.FacetGrid(data=data, col="vhgene", col_wrap=2, legend_out=False) # , height=1.5, aspect=.5)
    g.map_dataframe(sb.kdeplot, x='location',hue = 'antibody',shade=False).add_legend()
    plt.tight_layout()
    plt.show()

    '''
    sb.set(context='paper', style='ticks')
    g = sb.FacetGrid(data=df, col="vhgene", col_wrap=2, legend_out=False) # , height=1.5, aspect=.5)
    g.map_dataframe(sb.histplot, x='location',hue = 'antibody', bins=100, legend=True)
    g.add_legend()
    plt.tight_layout()
    plt.show()
    '''
    exit()

def figure_2b(): 
    
    filename1 = 'NeutSeqData_C002-215_cleaned_aligned_with_predictors.csv' 
    filename2 = 'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
    df = pd.read_csv('../../output/estimator/' + filename2 )


    ''' KDE cutoffs for VH3 '''
    dfs = df.loc[(df.coefficient.abs() > 1e-6)] 
    #dfs = df 
    dfs = dfs.loc[dfs['location'] <= 96]
    dfs = dfs.loc[dfs['chain'] == 'H']
    dfs['mut'] = dfs.aa + dfs['location'].astype(str)

    #vh3 =  ['F26', 'I27', 'A52', 'A56', 'T56', 'Y57', 'V97', 'A98','G98']

    #dfs = dfs.loc[dfs.mut.isin(vh3)] 
    print (dfs.mut) 
    locations = dfs['location'].unique()
    locations = pd.Series(locations).sort_values().values

    ''' beneficial mutations ''' 
    dfb = dfs.loc[dfs.coefficient < 0]
    #dfb = dfs
    dfb['coefficient_ab']  = dfs.coefficient.abs()

    print(dfb.mut.unique()) 

    print(dfb.head())

    exit()
    #size_order = dfb.coefficient.sort_values(ascending=True)


    ''' get a random sixteen antibodies for the poster '''
    #abs = list(dfb.antibody_id.values)
    #ablist = random.sample(abs, 14) + ['C105', 'REGN10986']
    #print(ablist)
    #dfbs= dfb.loc[dfb.antibody_id.isin(ablist)]
    dfbs = dfb 
    wt =['#636363','#e6550d']
    plt.figure(figsize=(2.75,5.5))
    sb.set(context='paper', style='white')
    #sb.scatterplot(data=dfb, x='mut', y = 'antibody_id', size='coefficient_ab', size_order= size_order, hue ='wild_type', palette=wt, legend='full')

    sb.scatterplot(data=dfbs, x='mut', y = 'antibody_id', size='coefficient_ab', hue ='wild_type', palette=wt, legend=False)
    plt.xticks(rotation =90)
    plt.tight_layout() 
    plt.savefig('../../output/figs/estimator-vh353.png', dpi=300)
    plt.show()    

    print(dfb.antibody_id.unique()) 

                    
#figure_4c()
figure_2b()
#figure_4e()

#figure_4b()
