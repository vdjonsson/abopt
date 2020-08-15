import pandas as pd
import numpy as np
import seaborn as sb
import math
from Levenshtein import distance as levenshtein_distance
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd
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

def plot_prediction(output_filepath, filename, df, y_col, legend_name, xlabel_name, log = True):
    plt.figure(figsize=[3.5,2.5], dpi=300)
    df.sort_values(by=[y_col], inplace=True)
    palette = ['#636363', '#3182bd', '#de2d26']
    if np.all(df[y_col].values > 0) and np.all(df[y_col+'_predicted'].values > 0):
        plt.semilogy(df.loc[:,df.columns == y_col].values, 'o', color = palette[0], alpha =0.5)
        plt.semilogy(df.loc[:,df.columns == y_col+'_predicted'].values, 'o', color = palette[1], alpha =0.5)
        plt.ylabel('log('+legend_name+')')
        #plt.legend(['log '+legend_name, 'log predicted '+legend_name])
    else:
        plt.plot(df.loc[:,df.columns == y_col].values.flatten(), 'o', color = palette[0], alpha = 0.5)
        plt.plot(df.loc[:,df.columns == y_col+'_predicted'].values.flatten(), 'o', color = palette[2], alpha =0.5)
        plt.ylabel(legend_name)
        #plt.legend([legend_name, 'predicted '+legend_name])
    plt.xticks([])
    plt.xlabel(xlabel_name)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_predictors_plot.png', dpi=300)
    plt.close()

def plot_estimator(output_filepath, filename, estimator_df, xlabel_name, ylabel_name, main_col):
    estimator = estimator_df.loc[estimator_df[main_col]!=0]
    if ylabel_name == 'Importances':
        cutoff = importance_cutoff(estimator[main_col].values)
        estimator = estimator.loc[estimator[main_col] > cutoff]
    if ylabel_name == 'Coefficients':
        neg, pos = coefficient_cutoff(estimator[main_col].values)
        estimator = estimator.loc[np.logical_or(estimator[main_col] < neg, estimator[main_col] > pos)]
    plt.figure(figsize=[9,3], dpi=300)
    plt.bar(range(len(estimator)), estimator[main_col].values, tick_label = estimator.index)
    plt.xticks(rotation=90)
    plt.xlabel(xlabel_name)
    plt.ylabel(ylabel_name)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_bar_'+ylabel_name+'.png', dpi=300)
    plt.close()

def create_levenshtein_map(df, id_col, heavy_col, light_col):
    l_distances = pd.DataFrame(data = np.array([np.nan]*(len(df)**2)).reshape((len(df), len(df))), columns = df[id_col].values, index = df[id_col].values)
    full_l_distances = pd.DataFrame(data = np.array([np.nan]*(len(df)**2)).reshape((len(df), len(df))), columns = df[id_col].values, index = df[id_col].values)
    
    col_combns = combinations(df[id_col].values, 2)
    for cols in list(col_combns):
        seq1_heavy, seq1_light = df.loc[df[id_col] == cols[0]][[heavy_col, light_col]].values.flatten()
        seq2_heavy, seq2_light = df.loc[df[id_col] == cols[1]][[heavy_col, light_col]].values.flatten()
        l_dist = levenshtein_distance(seq1_heavy, seq2_heavy) + levenshtein_distance(seq1_light, seq2_light)
        l_distances[cols[0]][cols[1]] = l_dist
        full_l_distances[cols[0]][cols[1]] = l_dist
        full_l_distances[cols[1]][cols[0]] = l_dist
    
    full_l_distances[np.isnan(full_l_distances)] = 0
    return full_l_distances, l_distances

def plot_l_heatmap(output_filepath, filename, l_distances):
    plt.figure(figsize=[10,7])
    sb.heatmap(l_distances, mask = l_distances.isnull(), xticklabels = True, yticklabels = True)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.savefig(output_filepath+filename+'_heatmap.png', dpi=300)
    plt.tight_layout()
    plt.close()

def create_linkage(l_distances):
    l_dist_array = ssd.squareform(l_distances)
    l_linkage = hierarchy.linkage(l_dist_array)
    return l_linkage
    
def plot_clustermap(output_filepath, filename, l_linkage, l_distances, df, id_col, y_col, patient_col, gene_col, y_name, patient_name, gene_name):
    cmap = cm.get_cmap('YlGnBu')
    y_colors = pd.Series([list(cmap(y/max(df[y_col].values))) for y in df[y_col].values], index = df[id_col], name = y_name)
    
    # to create legend:
    mappable = cm.ScalarMappable(norm = colors.Normalize(vmin=0, vmax=max(df[y_col].values)), cmap = 'YlGnBu')
    fig,ax = plt.subplots(1,1, figsize=(3,1.5))
    cb = fig.colorbar(mappable, ax=ax, aspect = 5)
    cb.set_label(label = y_name, size=6)
    cb.ax.tick_params(labelsize=4)
    cb.outline.set_visible(False)
    plt.savefig(output_filepath+filename+'_ic_legend.png', dpi=400)
    plt.close()
    
    if patient_col is not None:
        patient_colors_dict = dict(zip(np.unique(df[patient_col]), sb.color_palette('cubehelix', len(np.unique(df[patient_col])))))
        patient_colors = df[patient_col].map(patient_colors_dict)
        patient_colors.index = df[id_col]
        patient_colors.name = patient_name
        
        # to create legend:
        fig = plt.figure(figsize=(2, 3))
        patches = [
            mpatches.Patch(color=color, label=label)
            for label, color in patient_colors_dict.items()]
        fig.legend(patches, list(patient_colors_dict.keys()), loc='center', frameon=False, title = patient_name)
        plt.savefig(output_filepath+filename+'_patient_legend.png', dpi=300)
        plt.close()
        
    gene_colors_dict = dict(zip(np.unique(df[gene_col]), sb.color_palette('Paired', len(np.unique(df[gene_col])))))
    gene_colors = df[gene_col].map(gene_colors_dict)
    gene_colors.index = df[id_col]
    gene_colors.name = gene_name
    
    # to create legend:
    fig = plt.figure(figsize=(2, 1.25))
    patches = [
        mpatches.Patch(color=color, label=label)
        for label, color in gene_colors_dict.items()]
    fig.legend(patches, list(gene_colors_dict.keys()), loc='center', frameon=False, title = gene_name)
    plt.savefig(output_filepath+filename+'_gene_legend.png', dpi=300)
    plt.close()
    
    if patient_col is not None:
        g = sb.clustermap(l_distances, row_linkage = l_linkage, col_linkage = l_linkage, row_colors = pd.concat([y_colors, gene_colors], axis=1), col_colors = patient_colors, cbar_kws = {'label': 'edit distance'}, figsize=[6,6])
    else:
        g = sb.clustermap(l_distances, row_linkage = l_linkage, col_linkage = l_linkage, row_colors = y_colors, col_colors = gene_colors, cbar_kws = {'label': 'edit distance'}, figsize=[6,6])
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_clustermap.png', dpi=300)
    plt.close()

def plot_mapped_coefficients(filepath, location_filepath, output_filepath, filename, mapped_df, antibody_name, id_col):
    mapped_antibody = mapped_df.loc[mapped_df.index.get_level_values(id_col) == antibody_name]
    mapped_antibody = mapped_antibody.droplevel(id_col)
    neg, pos = coefficient_cutoff(mapped_antibody.coefficient.values)
    mapped_antibody = mapped_antibody.loc[np.logical_or(mapped_antibody.coefficient > pos, mapped_antibody.coefficient < neg)]
    
    locations = pd.read_csv(location_filepath+antibody_name+'_plotting_locations.csv', sep=',', header=0)
    pdb_locations = locations.pdb_location[list(mapped_antibody.index.get_level_values('location'))]
    mapped_antibody['pos'] = [aa+chain+str(pdb_location) for aa, chain, pdb_location in zip(list(mapped_antibody.index.get_level_values('aa')), list(mapped_antibody.index.get_level_values('chain')), pdb_locations)]
    
    g = sb.barplot(x='pos', y='coefficient', data = mapped_antibody, dodge=False, hue = 'wild_type', palette = {True: sb.color_palette('Set1')[0], False: sb.color_palette('Set1')[1]}, alpha = 0.7)
    g.legend().set_visible(False)
    plt.xlabel('Positions')
    plt.ylabel('Coefficients')
    plt.title(antibody_name)
    plt.xticks(rotation=90)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_'+antibody_name+'_mapped_coef_bar.png', dpi=300)
    plt.close()
    
def create_combined_df(dfs, id_cols, heavy_cols, light_cols, y_cols, patient_cols, gene_cols, new_id_col, new_heavy_col, new_light_col, new_y_col, new_patient_col, new_gene_col):
    new_df = pd.DataFrame(columns = [new_id_col, new_heavy_col, new_light_col, new_y_col, new_patient_col, new_gene_col])
    new_dfs = []
    for df, id_col, heavy_col, light_col, y_col, patient_col, gene_col in zip(dfs, id_cols, heavy_cols, light_cols, y_cols, patient_cols, gene_cols):
        if patient_col is None:
            df[new_patient_col] = 'NA'
            overlapping_ids = np.intersect1d(list(dfs[0][id_cols[0]]), list(df[id_col]))
            df = df[~df[id_col].isin(overlapping_ids)]
        else:
            df.rename(mapper={patient_col: new_patient_col}, axis=1, inplace=True)
        df.rename(mapper={id_col: new_id_col, heavy_col: new_heavy_col, light_col: new_light_col, y_col: new_y_col, gene_col: new_gene_col}, axis=1, inplace=True)
        df = df[[new_id_col, new_heavy_col, new_light_col, new_y_col, new_gene_col, new_patient_col]]
        new_dfs.append(df)
    new_df = pd.concat(new_dfs, join = 'inner', axis=0)
    return new_df
    

def plot_virus_wildtype(filepath, output_filepath, filename, estimator, true_sequence, offset):
    cutoff = importance_cutoff(estimator.importances.values)
    estimator = estimator.loc[estimator.importances > cutoff]
    estimator['wild_type'] = [False] * len(estimator)
    estimator['pos'] = estimator.index
    for i, position in enumerate(estimator.pos):
        aa = position[:1]
        index = int(position[1:])
        if true_sequence[index+offset] == aa:
            estimator.wild_type[i] = True
    plt.figure(figsize=(6.5,2.5))
    g = sb.barplot(x='pos', y='importances', data= estimator, dodge = False, hue = 'wild_type', palette = 'Set1', alpha=0.7)
    g.legend().set_visible(False)
    plt.xlabel('Positions')
    plt.ylabel('Importances')
    plt.title('SARS-CoV-2')
    plt.xticks(rotation=90)
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_wt_bar.png', dpi=300)
    plt.close()

def plot_logoplot(output_filepath, filename, coef_posmap, binding_sites, xlabel_name, ylabel_name, output_name, num_subplots, true_sequence):
    logomap = coef_posmap.T
    logomap.index = logomap.index.astype(int)
    fig, axs = plt.subplots(num_subplots, 1, figsize=(5.5,2.5))
    partition_size = math.floor(len(logomap)/num_subplots)
    for i in range(num_subplots):
        if i+1<num_subplots:
            sub_logomap = logomap[i*partition_size:(i+1)*partition_size]
            sub_true_sequence = true_sequence[i*partition_size:(i+1)*partition_size]
        else:
            sub_logomap = logomap[i*partition_size:]
            sub_true_sequence = true_sequence[i*partition_size:]
        logo = logomaker.Logo(sub_logomap, shade_below=0.5, fade_below=0.5, show_spines = False, ax = axs[i], color_scheme = '#636363')
        logo.style_glyphs_in_sequence(sequence=sub_true_sequence, color='#e6550d')
        for site in np.intersect1d(binding_sites, list(logomap.index[range(i*partition_size,(i+1)*partition_size)])):
            logo.highlight_position(p=int(site), color='#fff7bc', alpha = 0.5)
        if i+1 == num_subplots:
            logo.ax.set_xlabel(xlabel_name)
        logo.ax.set_ylabel(ylabel_name)
        logo.style_xticks(spacing=5)
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_'+output_name+'_logo.png', dpi=300)
    plt.close()
    
def plot_coef_heatmap(output_filepath, filename, coef_posmap):
    if np.any(coef_posmap > 0) and np.any(coef_posmap < 0):
        colormap = 'PuOr'
        center = 0
    else:
        colormap = 'PuBuGn'
        center = None
    f, ax = plt.subplots(figsize=(10, 5))
    ax = sb.heatmap(coef_posmap, cmap = colormap, center = center)
    plt.savefig(output_filepath+filename+'_coef_heatmap.png', dpi=300)
    plt.close()

def plot_hierarchical_clust(output_filepath, filename, l_distances, l_linkage, n_colors):
    hierarchy.set_link_color_palette([colors.rgb2hex(rgb[:3]) for rgb in cm.get_cmap('gist_earth')(np.linspace(0.1,0.9, n_colors))])
    dn = hierarchy.dendrogram(l_linkage, labels = l_distances.columns.values.tolist(), above_threshold_color = 'gray')
    plt.title('Hierarchical Clustering of Antibodies')
    plt.ylabel('Levenshtein distance')
    plt.tight_layout()
    sb.despine()
    plt.savefig(output_filepath+filename+'_dendrogram.png', dpi=300)
    plt.close()

def plot_violin_plot_mutations(output_filepath, filename, df, estimator):
    cutoff = importance_cutoff(estimator.importances.values)
    estimator = estimator.loc[estimator.importances > cutoff]
    mutations_estimator = pd.DataFrame(columns = ['residue', 'mutant', 'bind_avg'])
    for i,residue in enumerate(estimator.index):
        aa = residue[:1]
        pos = int(residue[1:])
        bind_avgs = df[['mutant', 'bind_avg']][np.logical_and(df.wildtype == aa, df.site_SARS2 == pos)]
        bind_avgs['residue'] = residue
        mutations_estimator = mutations_estimator.append(bind_avgs)
    # remove wild type
    mutations_estimator[mutations_estimator == 0] = np.nan
    plt.figure(figsize=(6.5,3))
    sb.violinplot(x = 'residue', y = 'bind_avg', data = mutations_estimator, color = 'white', scale = 'width', inner = None, linewidth = 0.5)
    g = sb.stripplot(x = 'residue', y = 'bind_avg', hue = 'mutant', data = mutations_estimator, palette= 'twilight', alpha=0.5,  dodge=False)
    # plt.legend(title = 'Amino acid', frameon = False, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    g.legend().set_visible(False)
    plt.xticks(rotation=90)
    plt.xlabel('Positions')
    plt.ylabel('$\Delta log(K_{d,ACE2})$')
    plt.tight_layout()
    plt.savefig(output_filepath+filename+'_violin.png', dpi=300)
    plt.close()


# nussenzweig paper

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../paper_figs/'
filename = 'NeutSeqData_C002-215_cleaned_aligned'
id_col = 'antibody_id'
heavy_col = 'igh_vdj_aa'
light_col = 'igl_vj_aa'
y_col = 'sars_cov_2_ic50_ngml'
patient_col = 'participant_id'
gene_col = 'heavy_v_gene'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
mapped_coefficients = pd.read_csv(filepath+filename+'_mapped_coefficients.csv', sep=',', header=0, index_col=[0,1,2,3])

plot_prediction(output_filepath, filename, metadata, y_col, '$IC_{50}(ng/ml)$', 'Antibody')
plot_estimator(output_filepath, filename, coefficients, 'Pseudopositions', 'Coefficients', 'coefficients')

coef_posmap = create_coef_matrix(coefficients)
plot_coef_heatmap(output_filepath, filename, coef_posmap)
plot_mapped_coefficients(filepath, location_filepath, output_filepath, filename, mapped_coefficients, 'C105', id_col)
specific_coefficients = mapped_coefficients.loc[np.logical_and(mapped_coefficients.index.get_level_values(id_col)=='C105', mapped_coefficients.index.get_level_values('chain') == 'H')]
wt_seq = ''.join(list(specific_coefficients[specific_coefficients['wild_type']==True].index.get_level_values('aa')))
specific_coefficients = specific_coefficients[['coefficient']]
specific_coefficients.index = [aa+str(location) for name,location,chain,aa in specific_coefficients.index]
specific_coef_map = create_coef_matrix(specific_coefficients)
binding_sites = pd.read_csv(location_filepath+'C105_contacts.csv', sep=',', header=0).number_antibody.values
plot_logoplot(output_filepath, filename, specific_coef_map, binding_sites, 'Positions', 'Coefficients', 'C105', 2, wt_seq)

full_l_distances, l_distances = create_levenshtein_map(metadata, id_col, heavy_col, light_col)
l_linkage = create_linkage(full_l_distances)
metadata[gene_col][~np.logical_or(metadata[gene_col].str.contains('HV3-53') , metadata[gene_col].str.contains('HV3-66'))] = 'Other'
plot_clustermap(output_filepath, filename, l_linkage, full_l_distances, metadata, id_col, y_col, patient_col, gene_col, '$IC_{50}$', 'patient', 'VH gene')
plot_hierarchical_clust(output_filepath, filename, full_l_distances, l_linkage, 9)

# VH3-53/66

combined_filename = filename + '_NeutSeqData_VH3-53_66_aligned'

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../paper_figs/'
filename = 'NeutSeqData_VH3-53_66_aligned'

vh3_metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
df = create_combined_df([metadata, vh3_metadata], [id_col, 'Name'], [heavy_col, 'VH or VHH'], [light_col, 'VL'], [y_col, 'IC50_ngml'], [patient_col, None], [gene_col, 'Heavy V Gene'], 'antibody_id', 'heavy_chain', 'light_chain', 'IC50', 'patient', 'VH gene')
full_l_distances, _ = create_levenshtein_map(df, 'antibody_id', 'heavy_chain', 'light_chain')
l_linkage = create_linkage(full_l_distances)
plot_clustermap(output_filepath, combined_filename, l_linkage, full_l_distances, df, 'antibody_id', 'IC50', 'patient', 'VH gene', '$IC_{50}$', 'patient', 'VH gene')

heavy_col = 'VH or VHH'
light_col = 'VL'
id_col = 'Name'
y_col = 'IC50_ngml'
patient_col = None
gene_col = 'Heavy V Gene'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)
mapped_coefficients = pd.read_csv(filepath+filename+'_mapped_coefficients.csv', sep=',', header=0, index_col=[0,1,2,3])

plot_prediction(output_filepath, filename, metadata, y_col, '$IC_{50}(ng/ml)$', 'Antibody')
plot_estimator(output_filepath, filename, coefficients, 'Pseudopositions', 'Coefficients', 'coefficients')

coef_posmap = create_coef_matrix(coefficients)
plot_coef_heatmap(output_filepath, filename, coef_posmap)
for antibody in ['C105', 'CB6', 'CV30', 'B38', 'CC12.1']:
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
    plot_logoplot(output_filepath, filename, specific_coef_map, binding_sites, 'Positions', 'Coefficients', antibody, 2, wt_seq)

full_l_distances, l_distances = create_levenshtein_map(metadata, id_col, heavy_col, light_col)
l_linkage = create_linkage(full_l_distances)
#plot_l_heatmap(output_filepath, filename, l_distances)
plot_clustermap(output_filepath, filename, l_linkage, full_l_distances, metadata, id_col, y_col, patient_col, gene_col, '$IC_{50}$', 'patient', 'VH gene')
plot_hierarchical_clust(output_filepath, filename, full_l_distances, l_linkage, 3)

# bloom paper

location_filepath = '../data/plotting/'
filepath = '../output/estimator/data/'
output_filepath = '../paper_figs/'
filename = 'single_mut_effects_cleaned'
y_col = 'bind_avg'
wt_seq = 'RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFHHHHHH'
binding_filename = 'SARS_CoV_2_ACE2_epitopes'

metadata = pd.read_csv(filepath+filename+'_with_predictors.csv', sep=',', header=0)
coefficients = pd.read_csv(filepath+filename+'_coefficients.csv', sep=',', header=0, index_col=0)

metadata[y_col] = -1*metadata[y_col]
metadata[y_col+'_predicted'] = -1*metadata[y_col+'_predicted']
plot_prediction(output_filepath, filename, metadata, y_col, r'$\Delta log(K_{d,ACE2})$', 'Virus mutant')
plot_estimator(output_filepath, filename, coefficients, 'Virus positions', 'Importances', 'importances')

coef_posmap = create_coef_matrix(coefficients)
plot_coef_heatmap(output_filepath, filename, coef_posmap)
plot_virus_wildtype(filepath, output_filepath, filename, coefficients, wt_seq, -1*(331-12))
ace_binding_sites = pd.read_csv(location_filepath+binding_filename+'.csv', sep=',', header=0).epitope_location.values
plot_logoplot(output_filepath, filename, coef_posmap, ace_binding_sites, 'Positions', 'Importances', 'sars_cov_2_ace2', 3, wt_seq)

c105_binding_sites = pd.read_csv(location_filepath+'C105_contacts.csv', sep=',', header=0).number_virus.values
plot_logoplot(output_filepath, filename, coef_posmap, c105_binding_sites, 'Positions', 'Importances', 'sars_cov_2_c105', 3, wt_seq)

plot_violin_plot_mutations(output_filepath, filename, metadata, coefficients)
