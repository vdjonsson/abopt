import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.cm
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd
from Levenshtein import distance as levenshtein_distance
from sklearn.manifold import MDS
from itertools import combinations

def get_sequence(database, name):
    sequence = database.loc[database.Name == name,:]['VH or VHH'].values[0]+database.loc[database.Name == name,:]['VL'].values[0]
    return sequence
    
filepath = '../data/'
db_filename = 'CoV-AbDab_270620'
data_filename = 'kyratsous_neutralization_data_ngml_units'

database = pd.read_csv(filepath+db_filename+'.csv', header=0, sep=',')
df = pd.read_csv(filepath+data_filename+'.csv', header=0, sep=',')
regn_ic50 = df.iloc[0,1:]

antibodies_dict = {}

for column in df.columns[1:]:
    antibodies_dict[column] = get_sequence(database, column)

data_filename = 'nussenzweig_antibody_data_cleaned_with_alignments'
df = pd.read_csv(filepath+data_filename+'.csv', sep=',', header=0)

for i in range(0,len(df)):
    antibodies_dict[df.antibody_id.values[i]] = df.igh_vdj_aa.values[i]+df.igl_vj_aa.values[i]

#heat map
l_distances = pd.DataFrame(data=np.array([np.nan]*(len(antibodies_dict)**2)).reshape((len(antibodies_dict),len(antibodies_dict))), columns = list(antibodies_dict.keys()), index = list(antibodies_dict.keys()))

col_combns = combinations(list(antibodies_dict.keys()), 2)
for cols in list(col_combns):
    seq1 = antibodies_dict[cols[0]]
    seq2 = antibodies_dict[cols[1]]
    l_distances[cols[0]][cols[1]] = levenshtein_distance(seq1, seq2)
    l_distances[cols[1]][cols[0]] = levenshtein_distance(seq2, seq1)
    
l_distances[np.isnan(l_distances)]=0
l_dist_array = ssd.squareform(l_distances)
l_linkage = hierarchy.linkage(l_dist_array)

nuss_ic50 = pd.DataFrame(data = df.sars_cov_2_ic50_ngml.values, index = df.antibody_id)
ic50_vals = pd.concat([regn_ic50, nuss_ic50]).astype(float)
cmap = matplotlib.cm.get_cmap('YlGnBu')
ic50_colors = pd.Series([list(cmap(ic50/max(ic50_vals.iloc[:,0].values))) for ic50 in ic50_vals.iloc[:,0].values], index = ic50_vals.index)

combined_output_filename = 'nussenzweig_kyratsous'

sb.clustermap(l_distances, row_linkage = l_linkage, col_linkage = l_linkage, row_colors = ic50_colors, col_colors = ic50_colors, figsize=(15,15))
plt.savefig('../figs/'+combined_output_filename+'_levenshtein_distance_clustermap.png')
plt.close()

dn = hierarchy.dendrogram(l_linkage, labels = l_distances.columns.values.tolist(), orientation = 'bottom')
print(dn['ivl'])
print(dn['color_list'])
plt.title('Hierarchical Clustering of Antibodies')
plt.ylabel('Levenshtein distance')
plt.savefig('../figs/'+combined_output_filename+'_levenshtein_distance_dendrogram.png')
plt.close()
exit(1)

plt.plot(range(0, len(ic50_vals)), np.log(ic50_vals.T[dn['ivl']].T), 'o')
plt.xticks(ticks = range(0, len(ic50_vals)), labels=dn['ivl'], rotation = 90, fontsize=4)
plt.tight_layout()
plt.savefig('../figs/'+combined_output_filename+'_levenshtein_distance_scatter_ic50.png')
plt.close()

mds = MDS(n_components = 2)
l_embedding = mds.fit(l_distances).embedding_
plt.plot(l_embedding[:,0], l_embedding[:,1], 'o')
plt.savefig('../figs/'+combined_output_filename+'_levenshtein_distance_MDS.png')
plt.close()

# plt.figure(figsize=[15,7], dpi=300)
# sb.heatmap(l_distances, mask = l_distances.isnull(), xticklabels=True, yticklabels=True)
# plt.xticks(fontsize=6)
# plt.yticks(fontsize=6)
# plt.savefig('../figs/kyratsous_nussenzweig_antibody_levenshtein_distance.png', dpi='figure')
# plt.close()

# plt.figure(dpi=300)
# plt.bar(range(0, len(ic50_vals)), ic50_vals.iloc[:,0].values, tick_label=ic50_vals.index)
# plt.xticks(fontsize=4, rotation = 90)
# plt.savefig('../figs/kyratsous_nussenzweig_ic50_vals.png')
# plt.close()
