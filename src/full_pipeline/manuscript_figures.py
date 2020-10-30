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
    
    ab_names = abdf.antibody.astype(str).values
    ab_class = abdf['class'].astype(str).values

    ab_class_dict = dict(zip(ab_names, ab_class))

    antibody_virus_landscape = pd.read_csv('../../output/merge/rbd_ab_fitness.csv')    

    ma.umap_antibody_distance(antibody_virus_landscape)


figure2c()
