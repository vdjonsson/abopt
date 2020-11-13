import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
import plotutils as pu 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.stats as stats 
from sklearn import preprocessing
import logomaker
import manuscript_analysis as ma

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


def figure2c():

    print('figure 2c') 

    ''' Read ab configuration file '''
    abdf = pd.read_table('../../manuscript/antibody_list.txt', sep=',') 
    
    antibody_virus_landscape = pd.read_csv('../../output/merge/rbd_ab_fitness.csv')    

    abt = antibody_virus_landscape.transpose() 
    abs = antibody_virus_landscape.columns[1:].values
    print(abs)

    abdf = abdf.loc[abdf.antibody.isin(abs)]
    ab_class = abdf['class'].values 
    ab_names = abdf['antibody'].values 
    ab_class_dict = dict(zip(ab_names, ab_class))

    classes = [ab_class_dict[x] for x in abs]

    ''' Graph embedding antibody dimensionality reduction '''
    embedding = ma.umap_antibody_distance(antibody_virus_landscape)
    print (embedding)

    embedding['ab'] = abs
    embedding['ab_class'] = classes

    print (embedding)


    ''' Plot ''' 
    sb.set(context='paper', style='white', font_scale=1.2)
    #f,ax = plt.subplots(figsize=(3.5,3.5))
    f,ax = plt.subplots()
    sb.set(context='paper')
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue= 'ab_class', s=80, legend=False)
    sb.despine()
    plt.gca().set_aspect('equal', 'datalim')
    for j, lab in enumerate(embedding.ab):
        ax.annotate(lab, (embedding.dim1[j] - 0.1, embedding.dim2[j]+ 0.1 ))                
    plt.title('UMAP')
    plt.tight_layout()
    plt.savefig('../../output/figs/' + 'UMAP_virus_scan_WT_ab.png', dpi = 300)
    plt.show()


    

def figure4(): 

    merge_dir = '../../output/merge/'
    df  = pd.read_csv(merge_dir + 'rbd_ace2_ab_fitness.csv')
    df = df.drop('ACE2', axis=1)

    ''' Analysis '''
    df['location'] = df.mut.str[1:-1].astype(int)
    df = df.loc[df['location'] <520]


    abcols = ['location','C105','B38','CC121','C002-S1', 'C002-S2','C119','C121-S1','C144','COVA2-39','C135','C110' ,'REGN10987','REGN10933','CV30', 'CB6']
    
    df = df[abcols]

    print(df.columns)
    
    dfs = df.iloc[:,1:-1]
    dfs = dfs.fillna(0)

    mask = np.zeros_like(dfs)
    print(mask)
    mask = dfs == 500

    ''' Normalize l2 '''
    normalizer = preprocessing.Normalizer()
    dfn = pd.DataFrame(normalizer.transform(dfs), columns=dfs.columns, index=df.location)
    df = dfn

    ''' Calculate mean ''' 
    grouped = df.groupby('location').mean().reset_index()
    grouped = grouped.set_index('location').transpose()

    mask = np.zeros_like(grouped)
    mask = grouped == 0  # mean equal to zero signifies that the structure is missing

    sb.set(context='paper')
    plt.figure(figsize=(13,2.5))
    g = sb.heatmap(data=grouped,mask=mask, cmap='RdBu_r')
    g.set_facecolor('#756bb1')
    plt.title('ddg(ab/RBD)<0 normalized')
    plt.tight_layout()
    plt.savefig(merge_dir + 'hm'+ '.png', dpi=300)
    plt.show()


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
    plt.savefig(merge_dir + 'hm_count'+ '.png', dpi=300)
    plt.show()
    
#figure2c()

figure4()
