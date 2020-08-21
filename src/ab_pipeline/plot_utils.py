import pandas as pd 
import seaborn as sb 
import matplotlib.pyplot as plt 
import antibody_pipeline_vdj as ab 
import numpy as np 
import utils as u

context = 'paper' 
style='ticks'
dpi = 300
out_fig ='../../output/figs/'
out_tmp ='../../output/tmp/'
est_path = '../../output/estimator/'
data_path = '../../data/'
ml_path = '../../data/ml_data/'


def plot_correlation(data, fignum):

    dcor = data.corr('pearson')

    u.write(dcor)
    figsize = (10,10)
    sb.clustermap (dcor, square=True, figsize=figsize, cmap='Reds', alpha=0.6)
    plt.title('Fitness correlation')
    plt.tight_layout()
    plt.savefig(out_fig + 'fitness_correlation_fig_' + str(fignum))
    plt.show()




def plot_ddg_distribution(pdb1, pdb2, title):

    ddg = ab.calculate_ddg_bind(pdb1 , pdb2, title)
    sb.set(context=context, style='ticks')
    sb.distplot( ddg.ddg.values, kde=False, hist=True)
    plt.semilogy()
    plt.ylabel('log10(ddg)')
    plt.title(title)
    plt.show()


def plot_ddg_heatmap_structures(pdb1, pdb2, title):

    ddg = ab.calculate_ddg_bind(pdb1 , pdb2, title)

    ddg = ddg[['mut','ddg']]
    ddg = ddg.set_index('mut')
    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    sb.heatmap(ddg)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_heatmap(ddg, title):

    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    plt.figure(figsize=(5,3))
    #sb.heatmap(ddg, square=False)
    sb.clustermap(ddg, figsize=(10,3), cmap='RdBu' )
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fig + 'hm.png', dpi=dpi)
    plt.show()


# input to this should be a ddg file with the location parsed out into another column 
def plot_ddg_stripplot(ddgsm,x, y, hue, filtername = None, filterval = '', title='ddg', subset=False, epitopes=None, epitopetype=None):

    #color_epitopes = dict(zip(colors, types))
    if filtername != None:
        print('non null filterval')
        ddgsm = ddgsm.loc[ddgsm[filtername] == filterval]
        print(ddgsm)
        print('printed')
    
    fig, ax = plt.subplots(figsize=(12,3))
    sb.violinplot(data=ddgsm,x=x,y=y,hue='ab', color='white',scale='width',inner=None, linewidth=0.5)
    g = sb.stripplot(data=ddgsm,x=x,y=y, hue='ab',s=4, palette='Set1', alpha=0.4, dodge=False)
    for epitope in epitopes:
        [plt.axvline(e, color ='#fee391', alpha = 0.6) for e in epitope]

    plt.xticks(rotation=90)
    ax.get_legend().remove()

    plt.tight_layout()
    plt.title(filterval)
    plt.tight_layout()
    plt.savefig(out_fig + title + '_' + filterval + '_strip.png', dpi=dpi)
    plt.show()


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



def scatter_plot_ddgs(ddgs):
    for m in muts_C105:
        sb.set(context=context, style=style)
        
        plt.figure(figsize=(2,2))
        sb.scatterplot(x=ddgs[wt].values, y=ddgs[m].values, color='Blue', alpha=0.2)
    
        xmin = ddgs[wt].min()
        ymin = ddgs[m].min()
        e= 0.3
        plt.hlines(y=0, xmin=xmin-e, xmax=1 , linestyles='dashed')
        plt.vlines(x=0, ymin=ymin-e, ymax=1, linestyles='dashed')

        plt.xlim([xmin-e,1])
        plt.ylim([ymin-e,1])

        plt.ylabel(m)
        plt.xlabel(wt)
        sb.despine()
        plt.tight_layout()
        plt.savefig(out_fig+'scatter_' + m + '.png', dpi=dpi)
        plt.show()


def dist_plot_ddgs(melted):
    # distribution neut, deleterious, gain
    sb.set(context=context, style=style)
    g = sb.FacetGrid(melted, row = 'bind_type',col="ab", hue = 'bind_type', height=2, margin_titles=True, sharex=False)
    g = g.map(sb.distplot, 'ddg', kde=False, bins=20)
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    plt.semilogy()
    plt.tight_layout()
    plt.savefig(out_fig + 'dist_ddg_v.png' , dpi=dpi)
    plt.show()

def plot_scatter(ab_name, y, yval):
    sb.set(context=context, style=style)
    plt.figure(figsize=(3,3))
    sb.scatterplot(x=range(len(y[yval].values)), y=y[yval].values, alpha=0.4)
    plt.ylabel(ab_name + ' ' + yval)
    plt.xlabel('Points')
    sb.despine()
    plt.tight_layout()
    plt.savefig(out_fig+ ab_name +'_' + yval + '.png', dpi=dpi)
    plt.show()


def plot_scatter_joint(ab_name, data, x, y, hue, kind = 'reg',leg=False):

    sb.set(context=context, style=style)
    plt.figure(figsize=(4,3))

    sb.scatterplot(data=data, x=x, y=y, hue=hue,alpha=0.4)

    #sb.jointplot(x=x, y=y, data=data, kind=kind, alpha=0.4, color='r')
    #plt.hlines(y=0, xmin=-1, xmax=1 , linestyles='dashed', lw=0.5)
    #plt.vlines(x=0, ymin=-1, ymax=1, linestyles='dashed', lw=0.5)

    if leg == False: 
        plt.legend([])

    plt.ylabel(y)
    plt.xlabel(x)
    sb.despine()
    plt.tight_layout()
    plt.savefig(out_fig+ 'scatter_' + ab_name +'_coeff_ddgs.png', dpi=dpi)
    plt.show()

def plot_antibody_viral_fitness():

    f = 'single_mut_effects_cleaned.csv'

    bloom = pd.read_csv(ml_path + f)
    df = bloom [['mutation', 'bind_avg']]

    df['bind_avg_bloom'] = df.bind_avg.values
    bloom, bloomval = normalize_data(df,colname='bind_avg_bloom')
    
    plot_scatter('bloom', bloom, bloomval)
    ab_name = 'B38'
    pdb_name = '7c01_Repair'
    mut = ''
    less_ab = '_less_ab_Repair'
    ab_name_graph = ab_name + '_' + mut

    df = ab.calculate_ddg_bind(ab_name,pdb_name, pdb_name + less_ab, scan_type='virus', mut_name='')
    ddgs = convert_to_simple_mutation(df)
    ddgs['ddg_' + ab_name_graph] = ddgs.ddg_.values
    ddgs, ddgsval = normalize_data(ddgs,colname='ddg_' +ab_name_graph)

    plot_scatter(ab_name, ddgs, ddgsval)
    merged = bloom.merge(ddgs, on='mutation')
    merged.to_csv(out_tmp + 'tmp.csv')

    plot_scatter_joint(ab_name_graph, merged, y=bloomval, x=ddgsval,hue='mutation', leg=False, kind='scatter')

 
def bar_plot_ddgs(data, ab_name, figsize):

    plt.figure(figsize=figsize)
    sb.set(context='paper', style='ticks')
    sb.barplot(data=data, x = 'mut', y='ddg_bind',palette='RdBu_r')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.title('ddgs mutations')
    plt.savefig(out_fig+ ab_name +'ddgs.png')
    plt.show()

