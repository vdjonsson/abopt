import pandas as pd
import numpy as np
import seaborn as sb
import math
from Levenshtein import distance as levenshtein_distance
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd
from scipy import stats, signal
from itertools import combinations
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from PIL import ImageColor
import logomaker
from matplotlib import cm, colors
import matplotlib.patches as mpatches
from plot_utils import coefficient_cutoff, importance_cutoff
from seqparser import create_coef_matrix
#sb.set(rc = {'font.size': 6})
#plt.rcParams.update({'font.size': 6})
#plt.rcParams.update({'mathtext.default': 'regular'})

def plot_estimator(output_filepath, filename, estimator_df, xlabel_name, ylabel_name, main_col):
    estimator = estimator_df.loc[estimator_df[main_col]!=0]
    if ylabel_name == 'Importances':
        cutoff = importance_cutoff(estimator[main_col].values)
        estimator = estimator.loc[estimator[main_col] > cutoff]
    if ylabel_name == 'Coefficients':
        neg, pos = coefficient_cutoff(estimator[main_col].values)
        estimator = estimator.loc[np.logical_or(estimator[main_col] < neg, estimator[main_col] > pos)]
    plt.figure(figsize=[8,2], dpi=300)
    plt.bar(range(len(estimator)), estimator[main_col].values, tick_label = estimator.index, color = '#636363', alpha = 0.7)
    plt.xticks(rotation=90)
    plt.xlabel(xlabel_name)
    plt.ylabel(ylabel_name)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_bar_'+ylabel_name+'.png', dpi=300)
    plt.close()

def plot_mapped_coefficients(filepath, location_filepath, output_filepath, filename, mapped_df, antibody_name, id_col):
    mapped_antibody = mapped_df.loc[mapped_df.index.get_level_values(id_col) == antibody_name]
    neg, pos = coefficient_cutoff(mapped_antibody.coefficient.values)
    mapped_antibody = mapped_antibody.loc[np.logical_or(mapped_antibody.coefficient > pos, mapped_antibody.coefficient < neg)]
    
    locations = pd.read_csv(location_filepath+antibody_name+'_plotting_locations.csv', sep=',', header=0)
    pdb_locations = locations.pdb_location[list(mapped_antibody.index.get_level_values('location'))]
    mapped_antibody['pos'] = [aa+chain+str(pdb_location) for aa, chain, pdb_location in zip(list(mapped_antibody.index.get_level_values('aa')), list(mapped_antibody.index.get_level_values('chain')), pdb_locations)]
    
    positive_mapped_antibody = mapped_antibody.loc[mapped_antibody['coefficient'] > 0]
    plt.figure(figsize=(6.5,2.5))
    g = sb.barplot(x='pos', y='coefficient', data = positive_mapped_antibody, dodge=False, hue = 'wild_type', palette = {True: sb.color_palette('Set1')[0], False: sb.color_palette('Set1')[1]}, alpha = 0.7)
    g.legend().set_visible(False)
    plt.xlabel('Positions')
    plt.ylabel('Coefficients')
    plt.title(antibody_name)
    plt.xticks(rotation=90)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_'+antibody_name+'_mapped_coef_bar_positive.png', dpi=300)
    plt.close()
    
    negative_mapped_antibody = mapped_antibody.loc[mapped_antibody['coefficient'] <= 0]
    plt.figure(figsize=(6.5,2.5))
    g = sb.barplot(x='pos', y='coefficient', data = negative_mapped_antibody, dodge=False, hue = 'wild_type', palette = {True: sb.color_palette('Set1')[0], False: sb.color_palette('Set1')[1]}, alpha = 0.7)
    g.legend().set_visible(False)
    plt.xlabel('Positions')
    plt.ylabel('Coefficients')
    plt.title(antibody_name)
    plt.xticks(rotation=90)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_'+antibody_name+'_mapped_coef_bar_negative.png', dpi=300)
    plt.close()

def plot_logoplot(output_filepath, location_filepath, filename, estimator_df, coef_posmap, binding_sites, xlabel_name, ylabel_name, output_name, num_subplots, true_sequence, main_col, mapped_df, antibody_name, id_col):
    logomap = coef_posmap.T
    estimator = estimator_df.loc[estimator_df[main_col]!=0]
    if ylabel_name == 'Importances':
        cutoff = importance_cutoff(estimator[main_col].values)
    elif ylabel_name == 'Coefficients':
        mapped_antibody = mapped_df.loc[mapped_df.index.get_level_values(id_col) == antibody_name]
        mapped_antibody = mapped_antibody.droplevel(id_col)
        neg, pos = coefficient_cutoff(mapped_antibody.coefficient.values)
    new_true_sequence = true_sequence
    logomap.index = logomap.index.astype(int)
    
    locations = pd.read_csv(location_filepath+antibody_name+'_plotting_locations.csv', sep=',', header=0)
    
    positive_logomap = logomap.clip(lower=0)
    positive_true_sequence = ''.join(np.array(list(new_true_sequence))[np.any(positive_logomap.values > pos, axis=1)])
    positive_logomap = positive_logomap.loc[np.any(positive_logomap.values > pos, axis=1)]
    indices = positive_logomap.index.values
    positive_logomap.index = list(range(len(positive_logomap.index.values)))
    fig, axs = plt.subplots(num_subplots, 1, figsize=(5.5,3.5))
    logo = logomaker.Logo(positive_logomap, shade_below=0.5, fade_below=0.5, show_spines = False, ax = axs, color_scheme = '#636363', flip_below=False)
    logo.style_glyphs_in_sequence(sequence=positive_true_sequence, color='#e6550d')
    logo.ax.set_xlabel(xlabel_name)
    logo.ax.set_ylabel(ylabel_name)
    tick_labels = [s+str(locations.pdb_location[i]) for s,i in zip(list(positive_true_sequence), indices)]
    logo.ax.set_xticks(positive_logomap.index.values)
    logo.ax.set_xticklabels(tick_labels, rotation=90)
    if len(binding_sites) > 0:
        logo.ax.set_ylim([np.min(positive_logomap.values), 1.1*np.max(positive_logomap.values)])
        logo.ax.set_yticks(np.linspace(np.min(positive_logomap.values), np.max(positive_logomap.values), 6).tolist())
        logo.style_spines(spines=['left'], visible=True, bounds=[np.min(positive_logomap.values), np.max(positive_logomap.values)])
        y = 1.1*np.max(positive_logomap.values)
        logo.ax.axhline(y, color='#FEBF5A', linewidth=10)
    for site in binding_sites:
        if len(np.where(indices == site)[0]) > 0:
            place = np.where(indices == site)[0].item()
            logo.ax.plot([place-0.5, place+0.5], [y,y], color = '#F43D25', linewidth=10, solid_capstyle='butt')
            
    plt.title(antibody_name)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_'+output_name+'_logo_positive.png', dpi=300)
    plt.close()
    
    negative_logomap = logomap.clip(upper=0)
    negative_true_sequence = ''.join(np.array(list(new_true_sequence))[np.any(negative_logomap.values < neg, axis=1)])
    negative_logomap = negative_logomap.loc[np.any(negative_logomap.values < neg, axis=1)]
    indices = negative_logomap.index.values
    negative_logomap.index = list(range(len(negative_logomap.index.values)))
    fig, axs = plt.subplots(num_subplots, 1, figsize=(5.5,3.5))
    logo = logomaker.Logo(negative_logomap, show_spines = False, ax = axs, color_scheme = '#636363', flip_below=False)
    logo.style_glyphs_in_sequence(sequence=negative_true_sequence, color='#e6550d')
    logo.ax.set_xlabel(xlabel_name)
    logo.ax.set_ylabel(ylabel_name)
    tick_labels = [s+str(locations.pdb_location[i]) for s,i in zip(negative_true_sequence, indices)]
    logo.ax.set_xticks(negative_logomap.index.values)
    logo.ax.set_xticklabels(tick_labels, rotation=90)
    
    if len(binding_sites) > 0:
        logo.ax.set_ylim([np.min(negative_logomap.values), np.max(negative_logomap.values) + 0.1*np.abs(np.min(negative_logomap.values))])
        logo.ax.set_yticks(np.linspace(np.min(negative_logomap.values), np.max(negative_logomap.values), 6).tolist())
        logo.style_spines(spines=['left'], visible=True, bounds=[np.min(negative_logomap.values), np.max(negative_logomap.values)])
        y = 0.1*np.abs(np.min(negative_logomap.values))
        logo.ax.axhline(y, color='#FEBF5A', linewidth=10)
    for site in binding_sites:
        if len(np.where(indices == site)[0]) > 0:
            place = np.where(indices == site)[0].item()
            logo.ax.plot([place-0.5, place+0.5], [y,y], color = '#F43D25', linewidth=10, solid_capstyle='butt')
            
    plt.title(antibody_name)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_'+output_name+'_logo_negative.png', dpi=300)
    plt.close()

# nussenzweig paper

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../new_paper_figs/'
filename = 'NeutSeqData_C002-215_cleaned_aligned'
id_col = 'antibody_id'
heavy_col = 'heavy_chain'
light_col = 'light_chain'
y_col = 'sars_cov_2_ic50_ngml'
patient_col = 'participant_id'
gene_col = 'heavy_v_gene'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
mapped_coefficients = pd.read_csv(filepath+filename+'_mapped_coefficients.csv', sep=',', header=0, index_col=[0,1,2,3])

coef_posmap = create_coef_matrix(coefficients)
for antibody in ['C105']: # c101, c135, c002
    plot_mapped_coefficients(filepath, location_filepath, output_filepath, filename, mapped_coefficients, antibody, id_col)
    specific_coefficients = mapped_coefficients.loc[np.logical_and(mapped_coefficients.index.get_level_values(id_col)==antibody, mapped_coefficients.index.get_level_values('chain') == 'H')]
    wt_seq = ''.join(list(specific_coefficients[specific_coefficients['wild_type']==True].index.get_level_values('aa')))
    specific_coefficients = specific_coefficients[['coefficient']]
    specific_coefficients.index = [aa+str(location) for name,location,chain,aa in specific_coefficients.index]
    specific_coef_map = create_coef_matrix(specific_coefficients)
    try:
        binding_sites = pd.read_csv(location_filepath+antibody+'_contacts.csv', sep=',', header=0).number_antibody.values
    except OSError as ose:
        binding_sites = []
    plot_logoplot(output_filepath, location_filepath, filename, coefficients, specific_coef_map, binding_sites, 'Positions', 'Coefficients', antibody, 1, wt_seq, 'coefficients', mapped_coefficients, antibody, id_col)

# VH3-53/66

combined_filename = filename + '_NeutSeqData_VH3-53_66_aligned'

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../new_paper_figs/'
filename = 'NeutSeqData_VH3-53_66_aligned'

heavy_col = 'heavy_chain'
light_col = 'light_chain'
id_col = 'antibody_id'
y_col = 'IC50_ngml'
patient_col = None
gene_col = 'Heavy V Gene'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
mapped_coefficients = pd.read_csv(filepath+filename+'_mapped_coefficients.csv', sep=',', header=0, index_col=[0,1,2,3])

coef_posmap = create_coef_matrix(coefficients)
for antibody in ['C105', 'CB6', 'CV30', 'B38', 'CC12.1']: # c101
    plot_mapped_coefficients(filepath, location_filepath, output_filepath, filename, mapped_coefficients, antibody, id_col)
    specific_coefficients = mapped_coefficients.loc[np.logical_and(mapped_coefficients.index.get_level_values(id_col)==antibody, mapped_coefficients.index.get_level_values('chain') == 'H')]
    wt_seq = ''.join(list(specific_coefficients[specific_coefficients['wild_type']==True].index.get_level_values('aa')))
    specific_coefficients = specific_coefficients[['coefficient']]
    specific_coefficients.index = [aa+str(location) for name,location,chain,aa in specific_coefficients.index]
    specific_coef_map = create_coef_matrix(specific_coefficients)
    try:
        binding_sites = pd.read_csv(location_filepath+antibody+'_contacts.csv', sep=',', header=0).number_antibody.values
    except OSError as ose:
        binding_sites = []
    plot_logoplot(output_filepath, location_filepath, filename, coefficients, specific_coef_map, binding_sites, 'Positions', 'Coefficients', antibody, 1, wt_seq, 'coefficients', mapped_coefficients, antibody, id_col)
# bloom paper

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../new_paper_figs/'
filename = 'single_mut_effects_cleaned'
y_col = 'bind_avg'
wt_seq = 'RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFHHHHHH'
binding_filename = 'SARS_CoV_2_ACE2_epitopes'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)

metadata[y_col] = -1*metadata[y_col]
metadata[y_col+'_predicted'] = -1*metadata[y_col+'_predicted']

coef_posmap = create_coef_matrix(coefficients)
ace_binding_sites = pd.read_csv(location_filepath+binding_filename+'.csv', sep=',', header=0).epitope_location.values
#plot_logoplot(output_filepath, location_filepath, filename, coefficients, coef_posmap, ace_binding_sites, 'Positions', 'Importances', 'sars_cov_2_ace2', 1, wt_seq, 'importances', mapped_coefficients, antibody, id_col)
