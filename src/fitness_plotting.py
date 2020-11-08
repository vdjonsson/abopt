import numpy as np
import pandas as pd
from sklearn.manifold import MDS, TSNE
import matplotlib.pyplot as plt

df = pd.read_csv('../data/rbd_ab_fitness.csv', sep=',', header=0, index_col=0)
df.dropna(inplace=True)
df = df.T

embedding = MDS(n_components=2).fit_transform(df.values)
fig, ax = plt.subplots()
ax.scatter(embedding[:,0], embedding[:,1])
for i, txt in enumerate(df.index.values):
	ax.annotate(txt, (embedding[i,0], embedding[i,1]))
plt.savefig('../paper_figs/MDS_rbd_ab_fitness.png', dpi=300)
plt.close()

embedding = TSNE(n_components=2).fit_transform(df.values)
fig, ax = plt.subplots()
ax.scatter(embedding[:,0], embedding[:,1])
for i, txt in enumerate(df.index.values):
        ax.annotate(txt, (embedding[i,0], embedding[i,1]))
plt.savefig('../paper_figs/TSNE_rbd_ab_fitness.png', dpi=300)
plt.close()
