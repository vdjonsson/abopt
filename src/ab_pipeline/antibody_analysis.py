import pandas as pd 
import seaborn as sb 
import matplotlib.pyplot as plt 
import antibody_pipeline_vdj as ab 


context = 'paper' 
#context = 'talk' 
style='ticks'
dpi = 300
out_fig ='../../output/figs/'
out_tmp ='../../output/tmp/'

def plot_ddg_distribution(p1, p2, title):

    ddg = ab.calculate_ddg_bind(p1 , p2, title)
    sb.set(context=context, style='ticks')
    sb.distplot( ddg.ddg.values, kde=False, hist=True)
    plt.semilogy()
    plt.ylabel('log10(ddg)')
    plt.title(title)
    plt.show()


def plot_ddg_heatmap_structures(p1, p2, title):

    ddg = ab.calculate_ddg_bind(p1 , p2, title)

    ddg = ddg[['mut','ddg']]
    ddg = ddg.set_index('mut')
    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    sb.heatmap(ddg)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_ddg_heatmap(ddg, title):

    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    plt.figure(figsize=(5,3))
    sb.heatmap(ddg, square=False)
    #sb.clustermap(ddg, figsize=(3,3))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fig + 'hm.png', dpi=dpi)
    plt.show()

def get_melted_ddgs(ddgs):
    ddgs = ddgs.reset_index()
    ddgs['location'] = ddgs.mut.str[3:-1]
    ddgsm = pd.melt(ddgs, id_vars=['mut', 'location'], value_vars=ddgs.columns[1:-1], var_name = 'ab',value_name='ddg')
    return ddgsm

def plot_ddg_stripplot(ddgsm,title='ddg', subset=False):


    plt.figure(figsize=(12,6))
    #sb.boxplot(data=ddgsm,x='location',y='ddg',palette='Paired')
    sb.violinplot(data=ddgsm,x='location',y='ddg',color='white',scale='width',inner=None, linewidth=0.5)
    sb.stripplot(data=ddgsm,x='location',y='ddg', hue='ab',s=5, palette='Paired', alpha=0.5, dodge=False)
    
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_fig + 'strip.png', dpi=dpi)
    plt.show()
    print(ddgsm.head())


def combine_data_sets(ddg_array):

    merged = pd.DataFrame()
    merged = ddg_array[0]

    for ddg_data in ddg_array[1:]: 
        print(ddg_data.head())
        merged = merged.merge(ddg_data, on='mut') 
    
    merged = merged.drop_duplicates(subset=['mut'], keep='first')
    merged = merged.set_index('mut')

    merged.to_csv(out_tmp + 'tmp.csv')
    return merged 

def get_subset_ddg(melted, yval = 'ddg', cutoff= 20):
    print(melted.head())
    melted = melted.loc[melted[yval]>=0]
    melted.to_csv(out_tmp + 'tmp.csv')
    y= melted.sort_values(by='ddg', ascending=False)
    return y 

def get_cutoff_ddg(melted, yval = 'ddg'):
    cutoff = 3 
    m = melted.loc[melted[yval]>=cutoff]
    y= m.sort_values(by='ddg', ascending=False)
    return y

def plot_ddg_distribution(melted, yval='ddg'):
    
    y = get_subset_ddg(melted, yval = 'ddg', cutoff= 0)

    cutoff = get_cutoff_ddg(melted)

    sb.set(context=context, style=style)
    plt.figure(figsize=(3,3))
    sb.scatterplot(x=range(len(y[yval].values)), y=y[yval].values, color='Blue', alpha=0.2)
    sb.scatterplot(x=range(len(cutoff[yval].values)), y=cutoff[yval].values, color='Red', alpha=0.2)
    
    plt.ylabel('ddg')
    plt.xlabel('points')
    sb.despine()
    plt.tight_layout()
    plt.savefig(out_fig+'dist.png', dpi=dpi)
    plt.show()
    

def bar_mutation_types(grouped):
    hue_colors = ['#6baed6','#9e9ac8', '#fb6a4a']
    hue_order = ['neg','neut','pos']

    # number neut, delet, gain mutations for these antibodies                                                       
    sb.set(context=context, style=style)
    sb.barplot(data=grouped, x='ab',y='ddg', hue='bind_type', hue_order = hue_order,palette=hue_colors)
    plt.xticks(rotation=90)
    plt.semilogy()
    plt.tight_layout()
    plt.savefig(out_fig + 'num_ddg_v_stacked.png' , dpi=dpi)
    plt.show()


# mutations of RBD in context of C105_wt

# mutations of RBD in context of C105_mt  
p1 = '6XCM_Repair_TH28D_Repair'
p2 = '6XCM_Repair_TH28D_Repair_less_ab_Repair' #less ab, mutate virus 

#plot_ddg_distribution(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')
#plot_ddg_heatmap(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')

ab_name = 'C105'
pdb_name = '6XCM_Repair'
less_str = '_less_ab_Repair'

mut_filename = ab_name + '_ab_mutations.txt'
muts = pd.read_table(ab.pdb_path + ab_name + '/' + mut_filename, header=None)
muts = muts.iloc[:,0].str.rstrip(';')


ddg_array = []
for mut in muts:

    mut_pdb = pdb_name + '_' + mut + '_Repair'
    print (mut_pdb)
    mut_name = pdb_name + '_' + mut
    ddg = ab.calculate_ddg_bind(ab_name, mut_pdb , mut_pdb + less_str, mut_name = mut_name) 
    ddg_array.append(ddg)

    
ab_name = 'B38'
pdb_name = '7bz5_Repair'
less_str = '_less_ab_Repair'

mut_pdb = pdb_name
mut_name = pdb_name
ddg = ab.calculate_ddg_bind(ab_name, mut_pdb , mut_pdb + less_str, mut_name = mut_name) 
ddg_array.append(ddg)


ab_name = 'CB6'
pdb_name = '7c01_Repair'
less_str = '_less_ab_Repair'

mut_pdb = pdb_name
mut_name = pdb_name
ddg = ab.calculate_ddg_bind(ab_name, mut_pdb , mut_pdb + less_str, mut_name = mut_name)
ddg_array.append(ddg)

ddgs = combine_data_sets(ddg_array)
ddgsm = get_melted_ddgs(ddgs)

#plot_ddg_heatmap(ddgs, title='ddg_' + ab_name)
print (ddgs.head())
#plot_ddg_stripplot(ddgsm,title='ddg', subset=False)

melted = get_melted_ddgs(ddgs)
#plot_ddg_distribution(melted, yval='ddg')


# just cutoff most postive
cutoff = get_cutoff_ddg(melted, yval = 'ddg')
#plot_ddg_stripplot(cutoff,title='ddg', subset=False)

# now look at antibody sensitivity to virus proportion of positive vs all for all designs 
print (ddgsm.head())

# determine mutation types: should be less stringent
melted.loc[melted.ddg > 0, 'bind_type'] = 'pos'
melted.loc[melted.ddg < 0, 'bind_type'] = 'neg'
melted.loc[(melted.ddg < 0.4) & (melted.ddg >-0.4), 'bind_type'] = 'neut'

all_muts = len(ddgsm.mut.unique())
grouped = melted.groupby(by=['bind_type','ab']).count()
grouped = grouped.reset_index()

hue_colors = ['#6baed6','#9e9ac8', '#fb6a4a']
hue_order = ['neg','neut','pos']

# bar_mutation_types(grouped)

# normalize pos and neg vals separately, cutoff values then concat structures  
# then heatmap 

maxval = ddgs.max().max()
minval = ddgs.min().min()


ddgs = ddgs.mask(ddgs > 0, ddgs/maxval)
ddgs = ddgs.mask(ddgs < 0, -ddgs/minval)
#ddgs = ddgs/minval


print(ddgs.head())
ddgs['sumval'] = ddgs.sum(axis=1)

sums = ddgs.sort_values('sumval', ascending=False)

cutoff = 0.4
cutoff_low = 0.1
sb.scatterplot(x=range(len(sums['sumval'].values)), y=sums['sumval'].values, color='Red', alpha=0.2)
plt.hlines(y=cutoff, color='b', xmin=0, xmax=len(sums['sumval'].values))
plt.hlines(y=-cutoff, color='b', xmin=0, xmax=len(sums['sumval'].values))
plt.ylabel('sum(ddg<0)')
plt.savefig(out_fig + 'sp_ddg_neg.png', dpi=dpi)
plt.show()

ddgs = ddgs.loc[(ddgs.sumval > cutoff) | (ddgs.sumval < -cutoff_low)]
print(ddgs.head())

ddgs = ddgs.drop(['sumval'], axis=1)
ddgs = ddgs.transpose()

sb.set(context=context, style=style)
plt.figure(figsize=(12,2))
sb.heatmap(ddgs, square=False,cmap='RdBu_r')
#sb.clustermap(ddg, figsize=(3,3), col_cluster =False)
plt.title('ddg(ab/RBD)<0 normalized')
plt.tight_layout()
plt.savefig(out_fig + 'hm_abs_neg.png', dpi=dpi)
plt.show()


exit()

# distribution neut, deleterious, gain
sb.set(context=context, style=style)
g = sb.FacetGrid(melted, row = 'bind_type',col="ab", hue = 'bind_type', height=2, margin_titles=True, sharex=False)
g = g.map(sb.distplot, 'ddg', kde=False, bins=20)
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
plt.semilogy()
plt.tight_layout()
plt.savefig(out_fig + 'dist_ddg_v.png' , dpi=dpi)
plt.show()


# 10933 6XDG
# now get subset 

