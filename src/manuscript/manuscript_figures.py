import sys
sys.path.insert(1, '/Users/vjonsson/Google Drive/data/repository/abopt-private/src/pipeline')

import colors as fmt 
import random 
import pandas as pd
import seaborn as sb 
import foldx as fx
import os 
import structure 
import energy as e     
import plotutils as pu 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import gcf
from sklearn.manifold import MDS, TSNE
import numpy as np 
import scipy.stats as stats 
import sklearn.preprocessing as pp
import logomaker
import manuscript_analysis as ma
import utils
from scipy import stats
import itertools as it 


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



'''  Figure 3 ''' 
def figure_3a(): 
    print ('Figure 3a')
    print ('Random forest estimator RBD')


def figure_3b(): 
    print ('Figure 3b')
    print ('Logo plot for random forest estimator')


def figure_3c(): 
    ''' Fitness landscape describing DDG(binding) of antibody  with respect to mutations of RBD 
        calculated by FoldX on molecular structures.  
    '''
    print ('Figure 3c')
    

def figure_3d():

    ''' UMAP of antibody RBD fitness with respect to point mutations of RBD   
    '''

    figname = 'figure_3d'

    alg = 'UMAP'

    ''' Read ab configuration file '''

    antibody_virus_landscape = utils.get_antibody_fitness_landscape(normalization=None)
    antibody_virus_landscape = antibody_virus_landscape.drop('C002-S2', axis=1) 

    abclass, abgene, ablabels =  utils.get_antibody_properties()

    abs = antibody_virus_landscape.columns[1:-1]
    abt = antibody_virus_landscape.iloc[:,1:-1]

    classes = [abclass[ab] for ab in abs]
    genes = [abgene[ab] for ab in abs]
    labels = [ablabels[ab] for ab in abs]

    ''' Graph embedding antibody dimensionality reduction '''
    embedding = ma.learn_antibody_distance(abt , alg = alg)
    embedding['ab'] = abs
    embedding['ab_class'] = classes
    embedding['ab_gene'] = genes
    embedding['ab_label'] = labels
    
    class_color = sb.color_palette('colorblind', 4)
    colors = class_color 
    hue = 'ab_class'

    area = np.random.rand(1, len(embedding))*40
    area = np.ones(len(embedding))*40

    ''' Plot ''' 
    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(3.75,3.75))
    sb.set(context='paper', font_scale=1)
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue, palette=colors, s=area, legend=False)
    plt.gca().set_aspect('equal', 'datalim')

    for j, lab in enumerate(embedding.ab_label):
        ax.annotate(lab, (embedding.dim1[j]-0.8, embedding.dim2[j]+ 0.3))                
    
    plt.title(alg)
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3d_' + alg + '_class.png', dpi = 300)
    plt.show()

    ''' Plot ''' 
    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(3,3))
    sb.set(context='paper', font_scale=0.7)
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue, palette=colors, s=area, legend=False)
    plt.gca().set_aspect('equal', 'datalim')

   # for j, lab in enumerate(embedding.ab_label):
    #    ax.annotate(lab, (embedding.dim1[j]-0.8, embedding.dim2[j]+ 0.3))                
    
    plt.title(alg)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3d_' + alg + '_class_nolabels.png', dpi = 300, transparent=True)
    plt.show()


def figure_2g_ab_scanning():
    
    ''' Fitness landscape for antibody design ''' 
    antibodies =['B38','COVA2-39','CC121'] #cc121 
    antibodies = ['CB6','CV30', 'C105']
    
    # ['C135', 'C144', 'C110', 'C119'] # more antibodies to add 

    ablist = pd.read_table('../../data/meta/antibody_list.txt', sep=',')
    ab_pdb = dict(zip(ablist.antibody.values, ablist.pdb.values))
    ab_hc = dict(zip(ablist.antibody.values, ablist.heavychain.values))
    
    filter = ['27','28','31','52','53','57','58', '76', '87','95','96']
    #filter = ['27']

    ddgall = pd.DataFrame()

    for ab in antibodies:
        repair_dir = '../../output/repair/'+  ab + '/' 
        pdb_file = ab_pdb[ab] +'_Repair' 
        pdb_file_less = ab_pdb[ab] +'_Repair_less_virus_Repair'
        scan_dir = '../../output/scan/'+ ab +'/'
        energy_dir = '../../output/energy/'+ ab +'/'
        location_file = '../../data/location/'+ ab + '_locations.csv'

        pos_scan_str = fx.construct_position_scan_string (pdb_name=ab_pdb[ab], location_file=location_file, chain = ab_hc[ab], filter= filter)
        fx.run_position_scan (pdb_file = pdb_file + '.pdb'  , scan_molecule = 'ab', pos_scan=pos_scan_str, pdb_dir = repair_dir, out_dir = scan_dir)    
        fx.run_position_scan (pdb_file = pdb_file_less + '.pdb' , scan_molecule = 'ab', pos_scan=pos_scan_str, pdb_dir = repair_dir, out_dir = scan_dir)
        ddg = e.calculate_ddg_bind(ab, pdb_file, pdb_file_less, scantype='ab', indir = scan_dir, outdir=energy_dir)
        ddgall = pd.concat([ddgall, ddg])

    ddgall.to_csv('../../output/merge/ddg_bind_ab_scanning.csv')


def figure_2g():


    figname = 'figure_2g_supp'

    antibodies =['B38',  'CV30', 'C105', 'CB6']


    ''' Merge the antibodies we want to graph ''' 

    ab_path = '../../output/energy/' 
    abs = [ab_path + ab + '/ddg_bind_' + ab + '_ab_scanning.csv' for ab in antibodies]
    
    abs_dict = dict(zip(antibodies, abs)) 

    data = pd.DataFrame()
    for ab in antibodies: 
        tmp = pd.read_csv(abs_dict[ab])
        data = pd.concat([tmp, data])

    data.to_csv('../../output/merge/ddg_bind_ab_scanning.csv')

    ''' Graph fitness landscape ''' 
    ablist = pd.read_table('../../data/meta/antibody_list.txt', sep=',')

    ab = pd.read_csv('../../output/merge/ddg_bind_ab_scanning.csv')
    ab['mutation'] = ab.mut.str[-1]
    ab = ab.loc[~ab.mutation.isin(['e','o'])]

    ''' Graph subset of antibody scanning data: not positions 52,76, and greater than 96 ''' 
    locs = [52, 76, 95,96]
    ab = ab.loc[~ab.pdb_location.isin(locs)]

    three = [ utils.aa_1to3_dic[x] for  x in ab.mutation.values]
    ab['mut_three'] = three


    data, design = pd.DataFrame(), pd.DataFrame()

    for antibody in antibodies: 
        abf = ab.loc[ab.antibody == antibody]        
        abf = abf[['mut_chain','pdb_location', 'mut_three', 'ddg_bind', 'antibody']] 

        ''' Calculate mean ddg and count ddg<0, for every pdb location ''' 
        design_tmp = abf.groupby('pdb_location').mean().reset_index()
        design_tmp['ddg_bind >=0']  = abf.loc[abf.ddg_bind >= 0].groupby('pdb_location').count().reset_index().ddg_bind.values 
        design_tmp['ddg_bind <0']  = 19*np.ones(len(design_tmp)) - design_tmp['ddg_bind >=0'].values 
        design_tmp['antibody'] = antibody 
        design_tmp['pdb_location'] = design_tmp.pdb_location.astype(str) 
        design = pd.concat([design, design_tmp])
    
        pivot = abf.pivot(index='mut_three', columns='pdb_location', values='ddg_bind')
        pivot = pivot.fillna(0)
        
        ''' Normalize l2 this should be changed '''
        normalizer = pp.Normalizer()
        data_norm = pd.DataFrame(normalizer.transform(pivot.values), columns=pivot.columns)

        data_norm ['mut_three'] = pivot.index.values
        data_norm['antibody'] = antibody
        #data_norm = data_norm.set_index('mut_three')

        data = pd.concat([data, data_norm])
    
    data = data.set_index('mut_three')
    area = design['ddg_bind <0'].values**1.2*4

    ''' Graph  mean ddg and count ddg<0, for every pdb location ''' 
    palette = 'Set1'
    plt.figure(figsize=(2.5,2.5))
    #plt.figure(figsize=(5,2.5))
    sb.set(context='paper', style='ticks', font_scale=1)
    sb.lineplot(data=design, x='pdb_location', y= 'ddg_bind', lw=0.6, ms=10, hue='antibody', legend=False, palette=palette) 
    sb.scatterplot(data=design, x='pdb_location', y= 'ddg_bind', lw=1.3,  hue='antibody', size=area, alpha=0.7, legend=True, palette=palette)
    plt.axhline(y=0, lw=0.5, color='grey') 
    plt.tight_layout()
    #plt.savefig(output_path + 'figure_2g.png', dpi=300, transparency=True)
    plt.savefig(output_path + 'figure_2g_legend.png', dpi=300, transparent=True)
    plt.show()
    exit()

    figsize = (3.5,4.5)
    figsize = (3.5,10)


    abs = data.pop("antibody")

    lut = dict(zip(abs.unique(), "rbgmk"))
    row_colors = abs.map(lut)

    sb.set(context='paper', style='ticks', font_scale=0.7)
    g = sb.clustermap(data=data,cmap='RdBu_r', center = 0,vmin=-1, vmax=1, row_cluster=False, col_cluster= True, figsize=figsize, row_colors=row_colors)

    print(abs)
    for label in abs.unique():
        print(label)
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],label=label, linewidth=0)

    g.ax_col_dendrogram.legend(loc="center", ncol=5)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_path + 'figure_2g_supp.png', dpi=300, transparent=True)
    plt.show()


    ''' Graph number of greater than zero and less than zero ddgs for each location '''
    

    exit()

def figure_4b():

    rbd = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_with_predictors.csv')

    ''' Read RBD fitness landscape ''' 
    rbd = rbd[['site_SARS2', 'bind_avg']]
    rbd['-bind_avg'] = -1 * rbd.bind_avg.values  
    vmin  = rbd['-bind_avg'].min()
    vmax  = rbd['-bind_avg'].max()

    ''' Read random forest features '''
    rf  = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_coefficients.csv')
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
    plt.savefig(output_path + 'figure_4b_ace2'+ '.png', dpi=300, transparent=True)
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
    plt.savefig(output_path + 'figure_4b'+ '.png', dpi=300, transparent=True)
    plt.show()

    ''' Plot '''
    sb.set(context='paper', style='ticks') 
    sb.violinplot(data=rbd, y='-bind_avg', x='site_SARS2', color='Red')
    plt.xticks(rotation=90)
    plt.show()


def calculate_C105_optimized_landscape(data):

    ''' Keep ddG binding for C105opt antibody only at interface of mutations and RBD
    '''

    ''' Find C105 epitope  '''
    epi = '6xcm_Repair_TH28I_YH58F_epitopes.csv'
    epdf  = pd.read_csv('../../output/epitope/'+ epi)

    epitope = epdf['number_virus'].unique()
    muts = data.loc[data['location'].isin(epitope)].mut.values
    
    data['C105_TH28I_YH58F_new'] = data.C105.values

    for mut in muts: 
        new_val = data.loc[data.mut==mut].C105_TH28I_YH58F.values[0]
        old_val = data.loc[data.mut==mut].C105.values[0]
        data.loc[data['mut'] == mut, 'C105_TH28I_YH58F_new'] = new_val

    data = data.drop('C105_TH28I_YH58F', axis=1)
    data = data.rename(columns={'C105_TH28I_YH58F_new':'C105_TH28I_YH58F'})

    return data

 


def initialize():

    merge_dir = '../../output/merge/'
    data = pd.read_csv(merge_dir + 'rbd_ab_fitness.csv')
    data_clean = clean_filter_landscape(data) 
    data_opt = calculate_C105_optimized_landscape(data_clean)
    
    data_opt.to_csv(merge_dir + 'rbd_ab_fitness_opt.csv', index=None)


def figure_3c(): 

    normalization = None 
    normalization_str = str(normalization)
    df = utils.get_antibody_fitness_landscape(normalization=normalization)

    ''' Calculate mean of ddgs at each PDB location, then add ab properties ''' 
    grouped = df.groupby('location').mean().reset_index()
    grouped = grouped.set_index('location').transpose()

    abclass, abgene, ablabels =  utils.get_antibody_properties()
    
    classv = [abclass[ab] for ab in grouped.index]
    genev = [abgene[ab] for ab in grouped.index]
    labelv = [ablabels[ab] for ab in grouped.index]

    grouped ['vh']  = genev
    grouped ['abclass'] = classv
    grouped ['label'] = labelv
    
    grouped = grouped.set_index('label')

    mask = np.zeros_like(grouped)
    mask = grouped == 0  # mean equal to zero signifies that the structure is missing

    ''' Add vhgene and antibody class ''' 
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


    ''' Get epitope locations '''     
    calculate_epitope = False
    if calculate_epitope:
        dir = '../../data/pdb/RBD_ACE2/'
        fname = '6m0j.pdb'
        distance = 6.0

        labeled_chains = structure.label_chains(fname[0:-4]) 
        epitopes = structure.find_epitopes(dir, fname, labeled_chains, distance)

    eps = pd.read_csv('../../data/location/epitopes.csv').epitope_location.values
    eps_rbd = []
    for loc in grouped.columns: 
        isep = False 
        if loc in eps: 
            isep = True 
        eps_rbd.append(isep)

    eps_labels = pd.Series(eps_rbd,index=grouped.columns)

    eps_pal = sb.color_palette('YlOrRd', eps_labels.unique().size)
    eps_lut = dict(zip(eps_labels.unique(), eps_pal))
    eps_colors = eps_labels.map(eps_lut)
    
    #col_colors = pd.DataFrame()
    col_colors = eps_colors
    col_colors = pd.DataFrame()
    col_colors['ACE2 epitope'] = eps_colors 

    grouped.to_csv('../../output/tmp.csv') 

    sb.set(context='paper', font_scale=0.65, style='ticks')
    figsize= (8,2.5) 
    fig, ax =plt.subplots()
    g = sb.clustermap(data=grouped, cmap='RdBu_r', vmin = -1, vmax = 1, center = 0, row_cluster = True, col_cluster=False,figsize=figsize, row_colors=row_colors, col_colors=col_colors, mask=mask , linecolor='white', linewidth=0.2,facecolor ='#bdbdbd')


    g.set_facecolor('#bdbdbd')
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
    plt.savefig(output_path + 'figure_4d_'+ normalization_str + '.png', dpi=300, transparent=True)
    plt.show()


    ''' Calculate number of positive ddgs '''
    dfs = dfs > 0
    #dfs['mut'] = df.mut
    dfs['location'] = df.index


    grouped = dfs.groupby(by='location').sum().reset_index()
    grouped = grouped.set_index('location').transpose()

    sb.set(context='paper')
    plt.figure(figsize=(12,2.5))
    g = sb.heatmap(data=grouped, mask=mask)#, cmap ='RdBu_r')
    g.set_facecolor('#756bb1')
    plt.tight_layout()
    plt.title('ddg(ab/RBD)<0 normalized')
    plt.tight_layout()
    plt.savefig(output_path + 'hm_count'+ '.png', dpi=300)
    #plt.show()

def melt_landscape(data):

    melted = data.melt(value_vars=data.columns[1:-1], id_vars=['location', 'mut'])
    melted = melted.rename(columns={'variable':'antibody', 'value':'binding'})
    return melted


def get_significant_binding_locations(antibody_1, antibody_2, test= 'ttest',sig_min=5e-2, fc_min=0):

    ''' Calculate differences in binding energies between any two antibodies  
        Returns dataframe with locations 
    '''
    abs = [antibody_1, antibody_2]
    data = utils.get_antibody_fitness_landscape(normalization='l2')
    melted = melt_landscape(data)
    filtered = melted.loc[melted.antibody.isin(abs)]
    classes, genes, labels = utils.get_antibody_properties()

    ''' Perform t test and compare accross locations designed ab vs wt  '''
    stat_test = pd.DataFrame()

    for loc in filtered['location'].unique():

        filtered_loc = filtered.loc[filtered['location'] == loc]
        ab1 = filtered_loc.loc[filtered_loc.antibody == abs[0]].binding
        ab2 = filtered_loc.loc[filtered_loc.antibody == abs[1]].binding

        ab1mean = ab1.mean()
        ab2mean = ab2.mean()

        res = stats.ttest_ind(ab1, ab2)          
        stat_test = stat_test.append(pd.Series([loc, res.pvalue, ab1mean, ab2mean]), ignore_index=True)

    ''' Extract statistically significant locations '''

    ab1_mean = abs[0] + '_mean' 
    ab2_mean = abs[1] + '_mean' 
    stat_test = stat_test.rename(columns={0:'location', 1:'pval', 2:abs[0] +'_mean', 3:abs[1] +'_mean'})
    stat_test['significant'] = stat_test.pval < sig_min
    stat_test['fold_change'] = np.abs(stat_test[ab1_mean]/stat_test[ab2_mean])

    sig = stat_test.loc[stat_test.significant] 

    print('Significant=' + str(len(sig)))
    sig = sig.loc[sig.fold_change > fc_min]
    print('Fold change > ' + str(fc_min) + '=' + str(len(sig)))

    sig_loc = sig.loc[sig.significant].location.astype(int).values

    filtered = filtered.loc[filtered['location'].isin(sig_loc)]
    filtered['significant'] = [x in sig_loc for x in filtered['location'].values]

    return filtered,sig

def get_significant_binding_locations_by_class(antibody_class_1, antibody_class_2, test= 'ttest',sig_min=5e-2, fc_min=0):

    ''' Calculate differences in binding energies between any two antibody classes   
        Returns dataframe with locations 
    '''

    abs1  = fmt.get_antibodies(antibody_class_1)
    abs2  = fmt.get_antibodies(antibody_class_2)

    data = utils.get_antibody_fitness_landscape(normalization=None)
    melted = melt_landscape(data)

    classes, genes, labels = utils.get_antibody_properties()

    filtered = melted 

    ''' Perform t test and compare accross locations designed ab vs wt  '''
    stat_test = pd.DataFrame()

    for loc in filtered['location'].unique():

        filtered_loc = filtered.loc[filtered['location'] == loc]

        ab_class1 = filtered_loc.loc[filtered_loc.antibody.isin(abs1)]
        ab_class2 = filtered_loc.loc[filtered_loc.antibody.isin(abs2)]

        ab1_class_mean= ab_class1.groupby('mut').mean().reset_index().binding.values
        ab2_class_mean = ab_class2.groupby('mut').mean().reset_index().binding.values

        ab1mean = ab1_class_mean.mean()
        ab2mean = ab2_class_mean.mean()

        res = stats.ttest_ind(ab1_class_mean, ab2_class_mean)          
        stat_test = stat_test.append(pd.Series([loc, res.pvalue, ab1mean, ab2mean]), ignore_index=True)


    ''' Extract statistically significant locations '''

    class_comp_str = antibody_class_1 + '_' + antibody_class_2
    ab1_mean_str = 'class_' + antibody_class_1 + '_mean' 
    ab2_mean_str = 'class_' + antibody_class_2 + '_mean' 
    stat_test = stat_test.rename(columns={0:'location', 1:'pval', 2:ab1_mean_str, 3:ab2_mean_str})
    stat_test['significant'] = stat_test.pval < sig_min
    stat_test['fold_change' ] = np.abs(stat_test[ab1_mean_str]/stat_test[ab2_mean_str])
    
    sig = stat_test.loc[stat_test.significant] 

    print('Significant=' + str(len(sig)))
    sig = sig.loc[sig['fold_change'] > fc_min]
    print('Fold change > ' + str(fc_min) + '=' + str(len(sig)))

    sig_loc = sig.loc[sig.significant].location.astype(int).values

    print(sig_loc) 

    filtered = filtered.loc[filtered['location'].isin(sig_loc)]
    filtered['significant'] = [x in sig_loc for x in filtered['location'].values]
    filtered['test_type'] = antibody_class_1 + '_' + antibody_class_2
 
    return filtered,sig


''' DEPRECATED '''
def figure_3e():
 
    figname = 'figure_3e'

    ''' Get list of combinations of two antibodies ''' 
    ablist = pd.read_table('../../data/meta/antibody_list.txt', sep=',')
    remove_abs = ['CR3022', 'C105_TH28I_YH58F', 'C104' , 'C121-S2', 'S309']
    ablist = ablist.loc[~ablist.antibody.isin(remove_abs)]
    abs = ablist.antibody.values 

    combinations = it.combinations(abs,2) 

    de = pd.DataFrame()
    for c in combinations: 

        print(c)
        abstr= c[0] + '_' +c[1]
        
        dfd = get_significant_binding_locations(antibody_1=c[0], antibody_2=c[1], test= 'ttest')

        tmp = pd.DataFrame()
        tmp['antibody_1'] =  c[0]
        tmp['antibody_2'] =  c[1]
        tmp['sig_location'] =  dfd.location.unique()

        de = pd.concat([de, tmp]) 

        ''' Plot locations '''
        plt.figure(figsize=(6.5,2.5))
        sb.set(context='paper', style='ticks', font_scale=1.3) 
        ax = sb.violinplot(data=dfd, x='location', y='binding',alpha=0.1,scale='width', hue='antibody',palette=fmt.c105opt, legend=True, split=True, inner='stick')
        ax.set_alpha(0.3)
        plt.axhline(0, lw=0.5, color = 'gray')
        plt.xticks(rotation=90)
        plt.title(c)
        ax.legend().set_visible(True)
        plt.tight_layout()   
        plt.savefig('../../output/manuscript/' + figname + '_' + abstr + '.png', dpi = 300, transparent=True)


    de.to_csv('../../output/merge/de_abd_rbd.csv')


def figure_3e_fitness():

    # All class 1 what are common mutations and what are specific to antibodies 
    class1 = all.loc[all['class'] =='1']
    class2 = all.loc[all['class'] =='2']
    class3 = all.loc[all['class'] =='3']

    class12 = class1.merge(class2, on='mutation', how='inner') 
    class13 = class1.merge(class3, on='mutation', how='inner') 
    class23 = class2.merge(class3, on='mutation', how='inner') 
    class123 = class12.merge(class3, on='mutation', how='inner')

    class123muts = class123.mutation.unique()
    class12muts = class12.mutation.unique()
    class13muts = class13.mutation.unique()
    class23muts = class23.mutation.unique()


    class1muts = class1.mutation.unique()

    classes = [class123muts, class12muts, class13muts, class23muts, class1muts]
    classes_name = ['classes123', 'classes12', 'classes13', 'classes23', 'class1']
    classes_list = [['1','2','3'], ['1','2'], ['1','3'], ['2','3'], ['1']]

    if False: 

        for i in range(len(classes)):

            c = classes[i]
            class_colors = fmt.get_class_colors(classes_list[i])
            print(class_colors) 
            smalldata = all.loc[all.mutation.isin(c)]
            smalldata = smalldata.loc[smalldata['class'].isin(classes_list[i])]

            sb.set(context='paper', style='ticks', font_scale=1)
            f,ax = plt.subplots(figsize=(7,2))
            plt.axhline(y=0, linewidth=0.5, color='gray')
            sb.stripplot(data = smalldata,x= 'location',y = 'AB_binding', hue='class', hue_order=classes_list[i],palette=class_colors, linewidth=0, alpha=0.6, marker='s')
            plt.gca().set_aspect('auto', 'datalim')
            plt.title('2D fitness landscape')
            plt.xlabel('ACE2_binding')
            plt.xticks(rotation=90)
            plt.ylabel('AB_binding')
            ax.legend().set_visible(False)
            plt.title(classes_name[i])
            plt.tight_layout()
            plt.savefig('../../output/manuscript/figure_3d_' + classes_name[i] + '.png', dpi = 300)
            plt.show()
            
            # Plot by antibody color 

            ab_list = smalldata.antibody.unique()
            ab_colors = fmt.get_antibody_colors(ab_list)
            
            sb.set(context='paper', style='ticks', font_scale=1)
            f,ax = plt.subplots(figsize=(7,2))
            plt.axhline(y=0, linewidth=0.5, color='gray')
            sb.stripplot(data = smalldata,x= 'location',y = 'AB_binding', hue='antibody', hue_order=ab_list,palette=ab_colors, linewidth=0, alpha=0.6, marker='s')
            
            plt.gca().set_aspect('auto', 'datalim')
            plt.title('2D fitness landscape')
            plt.xlabel('ACE2_binding')
            plt.xticks(rotation=90)
            plt.ylabel('AB_binding')
            ax.legend().set_visible(False)
            plt.title(classes_name[i])
            plt.tight_layout()
            plt.savefig('../../output/manuscript/figure_3d_antibody_' + classes_name[i] + '.png', dpi = 300)
            plt.show()


def figure_3c_new():

    figname = 'figure_c_new'

    f12, sig12 = get_significant_binding_locations_by_class('1', '2', test='ttest',sig_min=5e-2, fc_min=2)
    f13, sig13 = get_significant_binding_locations_by_class('1', '3', test='ttest',sig_min=5e-2, fc_min=2)
    f23, sig23 = get_significant_binding_locations_by_class('2', '3', test='ttest',sig_min=5e-2, fc_min=2)

    f12 = f12.loc[f12.significant].location.unique() 
    f13 = f13.loc[f13.significant].location.unique() 
    f23 = f23.loc[f23.significant].location.unique()

    alllocs = pd.Series(list(f12) + list(f13) + list(f23)).unique()

    print('=============')
    print(alllocs)
    print('=============')

    # Get fitness landscape, mean > 0 is escape potential C105, mean <=0 is C105 opt escape targeting 
    normalization = None
    normalization_str = str(normalization)

    data = utils.get_antibody_fitness_landscape(normalization=normalization)
    classes, genes, labels = utils.get_antibody_properties()
    melted = melt_landscape(data)

    classes = [classes[ab] for ab in melted.antibody.values]
    labels = [labels[ab] for ab in melted.antibody.values]
    melted['class'] = classes 

    ''' Get ACE2 binding '''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)

    # scale ace2 binding between -1 and 1 after adding + 0.1 minimum binding 
    print(ace2df) 


    ''' Merge ACE2 and AB binding ''' 
    merged = melted.loc[melted.location.isin(alllocs)]
    merged  = merged.merge(ace2df, how='inner', left_on='mut', right_on='mutation')

    cond = (merged.bind_avg >= fmt.MIN_ACE2_DDG_BIND) 
    merged = merged.loc[cond]  # escape mutants, better or same ACE2 binding and worse or same AB binding 
    merged['location']= merged.location.astype(str)

    grouped = merged.groupby(['location']).mean().reset_index()

    ''' Only graph antibodies with less than zero binding energies '''
    grouped = grouped.loc[grouped.binding.abs() >=0.01]
    locs = grouped.location

    merged = merged.loc[merged.location.isin(locs)]
    legend=False
    class_colors = fmt.get_class_colors(merged['class'].unique())
    hue_order = merged['class'].unique()
    merged['binding'] = merged.binding/merged.binding.max()

    print(merged.bind_avg.min())
    print(merged.bind_avg.max())


    ''' Mean binding '''
    hmdata = grouped[['location','bind_avg']].drop_duplicates().set_index('location').transpose()
    print(hmdata)
    plt.figure(figsize=(6,1.))
    sb.set(context='paper', font_scale=1., style='ticks')
    sb.heatmap(data=hmdata, cmap='Reds', linecolor='white', linewidth=0.3) 
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3e_new_muts_ic50_hm.png', dpi = 300, transparent=True)
    plt.show()

    ''' Plot mean binding of antibody class to mutations of RBD ''' 
    plt.figure(figsize=(6,2.5))
    sb.set(context='paper', font_scale=1., style='ticks')
    ax = sb.violinplot(data = merged, x= 'location', y = 'binding',color = 'white', alpha=0.4, jitter=1, dodge=False, scale= 'width', width =1., inner=None, linewidth=0.6)
    sb.stripplot(data = merged, x= 'location', y = 'binding',hue ='class', hue_order= hue_order, palette=class_colors,  alpha=0.8, jitter=1, dodge=False, size=3)
    plt.axhline(y=0, lw=0.5, color='gray')
    plt.title('2D fitness landscape')
    plt.xlabel('ACE2_binding')
    plt.xticks(rotation=90)
    plt.ylabel('AB_binding')
    #ax.legend().set_visible(legend)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3e_new_muts_mean.png', dpi = 300, transparent=True)
    plt.show()

    exit()

    ''' Now plot KDE's of several locations by class ''' 
    loc = 'all'
    locs = ['427', '440','484', '501' , '508']
    
    #locs = merged['location'].unique()

    variants = fmt.ALL_VARIANTS 

    exit()
    ''' This is not needed '''

    print(merged.columns)
    merged = merged[['location', 'antibody', 'class', 'binding', 'bind_avg', 'mut']]

    merged = merged.loc[merged.location.isin(locs)]

    merged.to_csv('../../output/manuscript/tmp.csv')
    class_colors = fmt.get_class_colors(merged['class'].unique())
    hue_order = merged['class'].unique()
        
    sb.set(context='paper', font_scale=1.2, style='ticks')
    grid = sb.FacetGrid(data=merged, col='class', row='location', hue='class',palette=class_colors, hue_order=hue_order, aspect=1)
    grid.map(sb.kdeplot,'bind_avg','binding', fill=True, alpha=0.8, legend=False)
    grid.map(sb.scatterplot,'bind_avg','binding',color='gray', marker='o', alpha=0.6)
    grid.map(plt.axhline,y=0, lw=0.5, color='gray')
    grid.map(plt.axvline,x=0, lw=0.5, color='gray')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname + '_' + str(loc) + '_kde_locs.png', dpi=300, transparent=True)
    plt.show()


    for location in locs: 
        mergedss = merged.loc[merged.location==location]

        print (mergedss)

        class_colors = fmt.get_class_colors(mergedss['class'].unique())
        hue_order = mergedss['class'].unique()
        
        sb.set(context='paper', font_scale=1.2, style='ticks')
        grid = sb.FacetGrid(data=mergedss, col='class', hue='class',palette=class_colors, hue_order=hue_order, aspect=1)
        grid.map(sb.kdeplot,'bind_avg','binding', fill=True, alpha=0.8, legend=False)
        grid.map(sb.scatterplot,'bind_avg','binding',color='gray', marker='o', alpha=0.6)
        grid.map(plt.axhline,y=0, lw=0.5, color='gray')
        grid.map(plt.axvline,x=0, lw=0.5, color='gray')
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_' + str(location) + '_kde_locs.png', dpi=300, transparent=True)
        plt.show()

    ''' Violin plots ''' 

    ''' heatmap ''' 

    pivot = merged.pivot_table(values='binding', index='mut', columns='antibody')
    plt.figure(figsize=(5,3))
    sb.set(context='paper', font_scale=1, style='ticks')
    ax = sb.clustermap(data = pivot, vmin=-1, vmax=1)
    plt.title('2D fitness landscape')
    plt.xlabel('ACE2_binding')
    plt.xticks(rotation=90)
    plt.ylabel('AB_binding')
    #ax.legend().set_visible(legend)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3e_new_muts.png', dpi = 300)
    plt.show()
 



def figure_3e():

    ''' Find locations on RBD with potential escape mutations for C105 and then color the impact of C105opt onto RBD protein.
    '''
    figname = 'figure_3e'
    annotate = True
    legend = False

    ''' Get all antibodies '''
    abs =['B38','CV30','CB6','C105','CC121', 'C119','C121','C144','COVA2-39','C135','C110' ,'REGN10987','REGN10933']

    # Get fitness landscape, mean > 0 is escape potential C105, mean <=0 is C105 opt escape targeting 
    normalization = None
    normalization_str = str(normalization)

    data = utils.get_antibody_fitness_landscape(normalization=normalization)
    classes, genes, labels = utils.get_antibody_properties()
    melted = melt_landscape(data)
    melted = melted.loc[melted.antibody.isin(abs)]

    classes = [classes[ab] for ab in melted.antibody.values]
    labels = [labels[ab] for ab in melted.antibody.values]

    melted['class'] = classes 

    ''' ACE2 binding '''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)

    all = pd.DataFrame()


    ''' Differential expression of binding to mutations per class of antibodies. 
        Take mean at each location per class, then compare with other classes.  
        Output list of locations that differ. 
    '''

    ''' Get all potential escape mutations that are infective ddG(ab) > 0 and binding ACE2 >=0 ''' 
    for ab in abs: 

        filtered = melted.loc[melted.antibody == ab]
        ab_escape  = filtered.loc[filtered.binding >0 ].sort_values('binding') 
        merged = ab_escape.merge(ace2df, left_on='mut', right_on ='mutation', how='inner') 

        cols = ['location', 'antibody','bind_avg', 'binding', 'mutation','class']
        renamed =  ['location', 'antibody','ACE2_binding', 'AB_binding', 'mutation','class']

        merged = merged[cols]
        merged = merged.rename(columns=dict(zip(cols, renamed)))
        merged = merged.loc[merged.ACE2_binding>=0]
        all = pd.concat([all, merged])

    hues = ['antibody', 'class']    

    ''' Plot every individual mutation ''' 
    for i in range(len(hues)):    

        hue = hues[i]

        ''' Group by location and calculate mean '''
        grouped  = all.groupby(['location' ,hue]).mean().reset_index()
        count = all.groupby(['location', hue]).count().reset_index()

        area = count.mutation.values**1*6

        ''' Get antibody and class colors ''' 
        if hue =='antibody':
            hue_color = fmt.get_antibody_colors(all.antibody.unique())
            hue_order = all.antibody.unique() 
        else:
            hue_color = fmt.get_class_colors(grouped['class'].unique())
            hue_order = grouped['class'].unique() 

        ''' Scale the binding data between 0 and 1 '''
        maxval = grouped.AB_binding.max()
        grouped['AB_binding'] = grouped.AB_binding/maxval

        maxval = grouped.ACE2_binding.max()
        grouped['ACE2_binding'] = grouped.ACE2_binding/maxval
        grouped = grouped.sort_values('location')
        grouped['location'] = grouped.location.astype(str)

        ''' Plot by location ''' 
        sb.set(context='paper', style='ticks', font_scale=1.2)
        f,ax = plt.subplots(figsize=(7,2))
        sb.scatterplot(data = grouped,x= 'location',y = 'AB_binding',hue = hue, hue_order=hue_order, palette=hue_color, legend=legend, alpha=0.8, s= area, marker='s')
        plt.gca().set_aspect('auto', 'datalim')
        plt.title('Potential escape variants, hue= ' + hue + ', size = number escape mutations')
        plt.xticks(rotation=90)
        plt.xlabel('PDB location')
        plt.ylabel('AB_binding')
        plt.tight_layout()
        plt.savefig('../../output/manuscript/figure_3e_grouped_location_' + hue +'_' + normalization_str + '.png', dpi = 300, transparent=True)
        plt.show()
        

        ''' Plot by binding ACE2 vs AB binding ''' 
        sb.set(context='paper', style='ticks', font_scale=1)
        f,ax = plt.subplots(figsize=(7,3.5))
        ax = sb.scatterplot(data = grouped,x= 'ACE2_binding',y = 'AB_binding',hue = hue,hue_order=hue_order,  palette=hue_color, legend=legend, alpha=0.8, s= area, marker='s')        
        plt.gca().set_aspect('auto', 'datalim')
        plt.title('2D fitness landscape')
        plt.xlabel('RBD location')
        plt.ylabel('AB_binding')
       
        c105opt= ['403', '405', '417', '420','421', '453','456','473' , '475', '476', '486', '495', '501','502','505']

        variants = fmt.ALL_VARIANTS

        if annotate:
            for j, lab in enumerate(grouped.location):
                bound1 = grouped.ACE2_binding[j] > 0.4
                bound2 = grouped.AB_binding[j] > 0.4
                bound3 = grouped.AB_binding[j] > 0.4
                bound4 = grouped.AB_binding[j] > 0.5
                bound5 = grouped.ACE2_binding[j] > 0.6
                cond = bound1 and bound2 or lab in c105opt
                #cond = (bound1 and bound2) or (bound3 and (lab in c105opt or lab in variants) or bound4 or bound5 and lab not in exclude)
                cond = lab in variants 

                if cond: 
                    ax.annotate(lab, (grouped.ACE2_binding[j]-0.008, grouped.AB_binding[j]+0.02))                

        plt.tight_layout()
        plt.savefig('../../output/manuscript/figure_3e_grouped_binding_' + hue + str(annotate) + '.png', dpi = 300, transparent=True)
        plt.show()

    ''' Plot 484 and also 501 '''
    mut1 = all.loc[all.mutation == 'E484K']
    mut2 = all.loc[all.mutation == 'N501Y']

    class_colors = fmt.get_class_colors(all['class'].unique())
    hue_order = all['class'].unique()
    all = all.sort_values('location')

    all['location'] = all['location'].astype(str)

    sb.set(context='paper', style='ticks', font_scale=1.2)
    f,ax = plt.subplots(figsize=(7,3))

    area = all.ACE2_binding**1*400

    all['AB_binding_scaled'] = all.AB_binding/all.AB_binding.max()
    print(area.min())
    print(area.max())
    print(all.ACE2_binding.min())
    print(all.ACE2_binding.max())

    
    sb.set(context='paper', font_scale=0.7)
    ax = sb.scatterplot(data = all, x= 'location', y = 'AB_binding',hue ='class', hue_order= hue_order, palette=class_colors,  alpha=0.9, s=area, legend=None)
    plt.title('2D fitness landscape')
    plt.xlabel('ACE2_binding')
    plt.xticks(rotation=90)
    plt.ylabel('AB_binding')
    ax.legend().set_visible(legend)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3d_all_muts' + normalization_str + '.png', dpi = 300)
    plt.show()

    # Get subsets of know escape mutations 

    exit()
    hue_colors = fmt.get_antibody_colors(abs)
    variants = fmt.CLASS2_ESCAPE_VARIANTS + fmt.CLASS3_ESCAPE_VARIANTS + fmt.CLASS12_ESCAPE_VARIANTS + fmt.CIRCULATING_VARIANTS 
    ss = all.loc[all.mutation.isin(variants)]

    hue_colors = fmt.get_antibody_colors(ss.antibody.unique())
    plt.figure(figsize=(4,4))
    sb.set(context='paper', style='ticks')
    ax = sb.scatterplot(data = ss, x= 'ACE2_binding', y = 'AB_binding_scaled',hue ='antibody', hue_order= hue_order, palette='Reds',  alpha=0.9, s=area, legend=True)
    #plt.gca().set_aspect('auto', 'datalim')    
    plt.title('2D fitness landscape')
    plt.xlabel('ACE2_binding')
    plt.xticks(rotation=90)
    plt.ylabel('AB_binding')
    ax.legend().set_visible(legend)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3d_all_muts' + normalization_str + '.png', dpi = 300, transparent=True)
    plt.show()



def figure_3f():

    ''' Plot significant binding energy differences between C105 and C105opt 
    with repect to mutations on RBD 
    '''
 
    figname = 'figure_3f'
    abs = ['C105_TH28I_YH58F','C105']

    filtered, sig = get_significant_binding_locations(antibody_1=abs[0], antibody_2=abs[1], test= 'ttest')
    print (filtered)

    ''' Restrict locations of epitope of where mutations where performed '''
    epi = '6xcm_Repair_TH28I_YH58F_epitopes.csv'
    epdf  = pd.read_csv('../../output/epitope/'+ epi)
    epitopes = epdf['number_virus'].unique()

    ''' Significant locations '''
    dfd = filtered.loc[filtered['location'].isin(epitopes)]
    legend = False

    abs = dfd.antibody.unique()
    ab_colors = fmt.get_antibody_colors(abs)
    hue_order = abs 

    ''' ACE2 binding '''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)

    merged = dfd.merge(ace2df, left_on='mut', right_on ='mutation', how='inner') 

    cond = (merged.bind_avg >= fmt.MIN_ACE2_DDG_BIND)
    merged_infective = merged[cond]

    grouped = merged_infective.groupby(['location']).mean().reset_index()

    print(grouped.bind_avg.min())
    print(grouped.bind_avg.max())

    ''' Mean binding '''
    hmdata = grouped[['location','bind_avg']].drop_duplicates().set_index('location').transpose()

    print(hmdata)
    plt.figure(figsize=(6,1.))
    sb.set(context='paper', font_scale=1., style='ticks')
    sb.heatmap(data=hmdata, cmap='Reds', linecolor='white', linewidth=0.3, vmin=-0.6, vmax=0.6)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3e_new_muts_ic50_hm.png', dpi = 300, transparent=True)
    plt.show()

    ''' Plot locations '''
    plt.figure(figsize=(6,3))
    sb.set(context='paper', style='ticks', font_scale=1.3) 
    ax = sb.violinplot(data=dfd, x='location', y='binding',scale='width', hue='antibody',palette=ab_colors, legend=False, split=True, inner='stick', lw=0.3, alpha=0.8, linewidth=0.5)
    ax.set_alpha(0.3)
    plt.axhline(0, lw=0.5, color = 'gray')
    plt.xticks(rotation=90)
    ax.legend().set_visible(legend)
    plt.tight_layout()   
    plt.savefig('../../output/manuscript/' + figname + '.png', dpi = 300, transparent=True)
    plt.show()



def figure_3g():

    ''' Find locations on RBD with potential escape mutations for C105 and then color the impact of C105opt onto RBD protein.
    '''
    figname = 'figure_3g'
    abs = ['C105_TH28I_YH58F','C105']

    # Get fitness landscape, mean > 0 is escape potential C105, mean <=0 is C105 opt escape targeting 
    data = utils.get_antibody_fitness_landscape(normalization='l2')
    melted = melt_landscape(data)
    filtered = melted.loc[melted.antibody.isin(abs)]

    filtered, sig = get_significant_binding_locations(antibody_1=abs[0], antibody_2=abs[1], test= 'ttest', fc_min=0)
    
    c105 = melted.loc[melted.antibody == abs[1]]
    c105 = c105.groupby('location').mean().reset_index()

    c105_escape = c105.loc[c105.binding >= 0.1 ].sort_values('binding')  #escape locations for C105 
    c105_escape['C105_escape'] = True 
    c105_c105opt = sig.sort_values('C105_TH28I_YH58F_mean') 
    
    ''' Merge the C105 escape map and the locations where binding of C105 and C105opt are significantly different ''' 
    merged = c105_escape.merge(c105_c105opt, how='left', on='location') 
    merged['C105opt_better'] = merged.C105_TH28I_YH58F_mean < 0 
    merged['C105opt_worse'] = merged.C105_TH28I_YH58F_mean >= 0 

    ''' Restrict locations of epitope of where mutations where performed '''
    epi = '6xcm_Repair_TH28I_YH58F_epitopes.csv'
    epdf  = pd.read_csv('../../output/epitope/'+ epi)
    epitopes = epdf['number_virus'].unique()

    isepitope = [e in epitopes for e in merged.location.values]
    merged['C105opt_epitope'] = isepitope 

    ''' Color ACE2 binding locations and random forest importances'''

    ace2 = 'SARS_CoV_2_ACE2_epitopes.csv'
    ace2df  = pd.read_csv('../../output/epitope/'+ace2)
    epitopes = ace2df['epitope_location'].unique()

    isepitope = [e in epitopes for e in merged.location.values]
    merged['ACE2_epitope'] = isepitope

    rf ='single_mut_effects_cleaned_coefficients.csv'
    rfdf  = pd.read_csv('../../output/estimator/'+rf)
    rfdf['location'] = rfdf.iloc[:,0].str[1:]

    grouped = rfdf.groupby('location').max().reset_index() 
    grouped['location'] = grouped.location.astype(int)
    print(rfdf.head())
    print(grouped.head())

    rfepi = grouped.loc[grouped.importances > 0.003]
    merged = merged.merge(rfepi, how='left', on='location')

    print(merged.head())
    merged.to_csv('../../output/map/C105_C105opt_resistance_map.csv') 


def figure_2c(): 

    ''' Figure 2c: dotplot of estimator coefficients < 0 
    '''
    
    filename1 = 'NeutSeqData_C002-215_cleaned_aligned_with_predictors.csv' 
    filename2 = 'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
    df = pd.read_csv('../../output/estimator/' + filename2 )

    print(df.head())
    

    ''' Use locations derived from cutoffs '''
    locations = [27, 28, 53, 57,58, 96,97]


    ''' KDE cutoffs for VH3 '''
    dfs = df.loc[(df.coefficient.abs() > 1e-6)] 
    dfs = dfs.loc[dfs['location'] <= 96]
    dfs['pdb_location'] = dfs['location'] + 1
    dfs = dfs.loc[dfs['chain'] == 'H']
    dfs['mut'] = dfs.aa + dfs['pdb_location'].astype(str)


    locations = dfs['location'].unique()
    locations = pd.Series(locations).sort_values().values

    ''' beneficial mutations ''' 
    dfb = dfs.loc[dfs.coefficient < 0]
    dfb['coefficient_ab']  = dfs.coefficient.abs()

    exclude = ['COV2-2113', 'COV2-2812','COV2-2752','COV2-2165','C102', 'COVA2-20']
    dfb = dfb.loc[~dfb.antibody_id.isin(exclude)]
    
    area = dfb.coefficient_ab**8*1e7

    plt.figure(figsize=(2.25,5.))
    sb.set(context='paper', style='ticks')
    ax =sb.scatterplot(data=dfb, x='mut', y = 'antibody_id', size=area, hue ='wild_type', palette=fmt.wildtype, legend=None , marker='s')
    plt.xticks(rotation = 90)
    plt.tight_layout() 
    plt.savefig(output_path +'estimator-vh353.png', dpi=300)
    plt.show()    

    ''' Legend '''
    figsize = (2.75, 2.75)
    fig_leg = plt.figure()
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.legend(*ax.get_legend_handles_labels(), loc='center')

    ax_leg.axis('off')
    fig_leg.savefig(output_path +'estimator-vh353-legend.png', dpi=300)
    plt.show()


def figure_3h():

    ''' Find locations on RBD with potential escape mutations for C105 and then color the impact of C105opt onto RBD protein.
    '''
    figname = 'figure_3h'

    annotate = True
    legend = True

    abs =['C105','C105_TH28I_YH58F']

    # Get fitness landscape, mean > 0 is escape potential C105, mean <=0 is C105 opt escape targeting 
    data = utils.get_antibody_fitness_landscape(normalization=None)
    cols =['location','mut'] + abs
    datass = data[cols]

    # significant locations 
    c105opt= ['403', '405', '417', '420','421', '453','456','473' , '475', '476', '486', '495', '501','502','505']
    variants = fmt.ALL_VARIANTS
    datass = datass.loc[datass['location'].isin(c105opt)]

    ''' ACE2 binding '''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)

    merged = datass.merge(ace2df, left_on='mut', right_on ='mutation', how='inner') 

    # Difference between optimized and non opt antibodies on binding of RBD mutations 
    merged['dddG'] = merged[abs[1]] - merged[abs[0]]
    merged['dddG_count<=0'] = merged.dddG <=0 

    # Frequency bar plot 
    grouped = merged.groupby('dddG_count<=0').count().reset_index()

    ab_colors = fmt.get_antibody_colors(abs)

    plt.figure(figsize=(2.,2.))
    sb.set(context='paper', font_scale=1, style='ticks')
    sb.barplot(data=grouped,x='dddG_count<=0',y='dddG', palette=ab_colors) 
    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname +'.png' , dpi=300)  
    plt.show()

    merged.to_csv('../../output/manuscript/tmp.csv')

    area = merged.bind_avg.abs()*500
    area = np.ones(len(merged.bind_avg))*50

    allloc = len(merged.location.unique())
    palette = sb.color_palette("Paired", 12) + sb.color_palette("Set2", allloc -12)  

    print (allloc)
    print (palette) 


    ''' Multiple kde plots ''' 

    variants = fmt.ALL_VARIANTS

    range_ab = [-1.5,1.5]
    range_ace2 = [-6,1.5]


    ssloc =['all','417','475','483','484','501'] 

    variants = fmt.ALL_VARIANTS
    colors = sb.color_palette('Set3', len(ssloc))

    i = 0

    for loc in ssloc: 

        mergedss = merged

        if loc != 'all':
            mergedss= merged.loc[merged.location == int(loc)]

        mergedvariants = mergedss.loc[mergedss.mut.isin(variants)]

        print(mergedvariants)
    

        plt.figure(figsize=(2.5,2.5))
        sb.set(context='paper', font_scale=1.2, style='ticks')
        sb.kdeplot(data=mergedss, x='bind_avg', y='dddG', fill=True, color=colors[i], alpha=0.8)
        sb.scatterplot(data=mergedss, x='bind_avg', y='dddG',color='gray', marker='o', alpha=0.6)
        ax =sb.scatterplot(data = mergedvariants,x= 'bind_avg', y='dddG',color='red', marker='o', alpha=0.6)
        plt.axhline( y=0, lw=0.5, color='gray')
        plt.axvline( x=0, lw=0.5, color='gray')
        plt.xlim(range_ace2)
        plt.ylim(range_ab)
        plt.title(str(loc))
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_' + loc + 'kde_locs.png', dpi=300, transparent=True)
        plt.show()
        i=i+1


    exit()


    legend = False
    plt.figure(figsize=(5,3))
    sb.set(context='paper', font_scale=1, style='ticks')
    ax =sb.scatterplot(data=merged, x='bind_avg',y = 'dddG', hue ='location',palette=palette,s=area, alpha=0.8)
    plt.axhline(y=0, lw=0.5, color='gray')
    plt.axvline(x=0, lw=0.5, color='gray')
    plt.xlabel('ACE2 binding')
    plt.ylabel('AB binding')

    if annotate:
        for j, lab in enumerate(merged.mut):
            cond = merged['mut'][j] in variants
            print (merged.mut[j])
            if cond:
                ax.annotate(lab, (merged['bind_avg'][j], merged['dddG'][j]))

    plt.xlim([-5,1])
    plt.xlabel('ddG bind (mutACE2)')
    plt.ylabel('dddG bind (ABopt/mutACE2, ABwt/mutACE2)')
    #plt.legend(ncol=5)
    ax.legend().set_visible(legend)
    plt.gca().set_aspect('auto', 'datalim')    
    plt.tight_layout()
    plt.savefig('../../output/manuscript/figure_3h_all_muts.png', dpi = 300)
    plt.show()
    plt.show()

    exit()


    area = merged.bind_avg.abs()**1.4*10
    ax =sb.scatterplot(data=merged, x=abs[1],y = abs[0], hue ='location',palette='Set1',s=area, alpha=0.8)
    #ax =sb.scatterplot(data=melted, x='ACE2_binding',y ='AB_binding', hue ='antibody',palette='Set1',s=area, alpha=0.8)
    plt.axhline(y=0, lw=0.5, color='gray')
    plt.axvline(x=0, lw=0.5, color='gray')
    if annotate:
        for j, lab in enumerate(merged.mut):
            cond = merged['mut'][j] in variants
            if cond:
                ax.annotate(lab, (merged[abs[1]][j], merged[abs[0]][j]))

    plt.xlim([-0.4,0.4])
    plt.gca().set_aspect('auto', 'datalim')    
    plt.savefig('../../output/manuscript/figure_3h_all_muts.png', dpi = 300)
    plt.show()
    plt.show()


    classes, genes, labels = utils.get_antibody_properties()
    melted = melt_landscape(data)
    melted = melted.loc[melted.antibody.isin(abs)]
    
        
    classes = [classes[ab] for ab in melted.antibody.values]
    labels = [labels[ab] for ab in melted.antibody.values]

    melted['class'] = classes 


    ''' ACE2 binding '''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)

    all = pd.DataFrame()

    ''' Get all potential escape mutations that are infective ddG(ab) > 0 and binding ACE2 >=0 ''' 
    for ab in abs: 

        filtered = melted.loc[melted.antibody == ab]
        #ab_escape  = filtered.loc[filtered.binding > 0].sort_values('binding') 
        ab_escape  = filtered

        merged = ab_escape.merge(ace2df, left_on='mut', right_on ='mutation', how='inner') 
        merged = merged.loc[merged['location'].isin(c105opt)]

        cols = ['location', 'antibody','bind_avg', 'binding', 'mutation','class']
        renamed =  ['location', 'antibody','ACE2_binding', 'AB_binding', 'mutation','class']

        merged = merged[cols]
        merged = merged.rename(columns=dict(zip(cols, renamed)))
        #merged = merged.loc[merged.ACE2_binding>=0]
        all = pd.concat([all, merged])

    print(all.head())

    hue_colors = fmt.get_antibody_colors(abs)
    hue_order = abs

    all = all.sort_values('location')
    all['location'] = all['location'].astype(str) 

    area = all.ACE2_binding.abs()**1.4*10
    plt.figure(figsize=(7,3))
    sb.set(context='paper', font_scale=1, style='ticks')
    plt.axhline(y=0, lw=0.5, color='gray')
    ax = sb.scatterplot(data = all, x= 'location', y = 'AB_binding',hue ='antibody', hue_order= hue_order, palette=hue_colors,  alpha=0.9, legend=True, s= area)
    plt.gca().set_aspect('auto', 'datalim')    
    plt.title('2D fitness landscape')
    plt.xlabel('ACE2_binding')
    plt.xticks(rotation=90)
    plt.ylabel('AB_binding')
    ax.legend().set_visible(legend)
    plt.tight_layout()

    c105opt= ['403', '405', '417', '420','421', '453','456','473' , '475', '476', '486', '495', '501','502','505']
    variants = ['E484K','Q493R', 'N439K','N440K', 'S477', 'V483A', 'N501Y']
    exclude = ['517', '471']

    grouped = all 

    grouped = grouped.loc[grouped.mutation.isin(variants)]

    print(grouped.head())
    exit()


''' Add robustness '''
def figure_5():

    allsims = '../../output/cocktail/cocktails_allsims_noise.csv'
    
    df = pd.read_csv(allsims) 

    df = df.groupby(['cov', 'noise']).min().reset_index()

    #covs = [0.5800000000000001, 0.81, 0.89, 0.91]
    #grouped = grouped.loc[grouped['cov'].isin(covs)]

    '''
    sb.violinplot(data=df, y="num_abs", x="cov", shade=True,palette='RdPu_r', legend=False)
    plt.show()

    plt.figure(figsize=(3,2))
    sb.kdeplot(data=df, x="num_abs", hue="cov", shade=True,palette='crest', legend=True,  alpha=.5, linewidth=0,)
    plt.xlabel('Number of antibodies')
    plt.tight_layout()
    plt.show()

    '''

    plt.figure(figsize=(3,2))
    
    #sb.lineplot(data=df, x="cov", y='num_abs', hue="noise", palette='crest', legend=True,  alpha=1., linewidth=1)
    sb.lineplot(data=df[df.noise=='noise_0'], x="cov", y='num_abs', hue="noise", color='Red', legend=True,  alpha=1., linewidth=1)
    sb.lineplot(data=df, x="cov", y='num_abs', hue='noise', palette='Reds', alpha=1, linewidth=0.5)
    plt.xlabel('Number of antibodies')
    plt.tight_layout()
    plt.show()

    ''' for all noise '''

    print(grouped.head())
    print(df.columns)
    


# set the output path

output_path  = '../../output/manuscript/'

# initialize C105 landscape first and save rbd_ace2 
#initialize() 
#figure_4c()
#figure_2c()
#figure_2c()
#figure_2g_ab_scanning()
#figure_3c()  
#figure_3d()  
#figure_3c_new()  
#figure_3f()  
#figure_3g()
#figure_3h()
#figure_4c()
#figure_4e()
#figure_4b()

figure_5()
