import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
from scipy import stats
import numpy as np 

from sklearn.manifold import MDS, TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt 
import umap 
from sklearn.model_selection import train_test_split

filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'

dpi = 300 
plot =  False
ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30', 'REGN10987']
ab_names_list = ['B38', 'C105','CC12.1', 'CB6', 'COVA2-39','CV30', 'REGN10987']


pdb_names = ['7bz5', '6xcm', '6xc3_CC121', '7c01', '7jmp', '6xe1' , '6xdg_REGN10987']
rbd_chains = ['A', 'C', 'C', 'A', 'A', 'E', 'E']
class_names= ['class I', 'class I', 'class I', 'class I', 'class II','class I', 'class III']


beg = 0
end = 7
ab_names = ab_names[beg:end]
pdb_names =pdb_names[beg:end]
class_names = class_names[beg:end]
rbd_chains = rbd_chains[beg:end]

#ab_names = [ab_names[1]]
#pdb_names =[pdb_names[1]]
#class_names = [class_names[1]]
#rbd_chains = [rbd_chains[1]]


'''
ab_names = [ab_names[1], ab_names[2],ab_names[4]]
pdb_names =[pdb_names[1],pdb_names[2],pdb_names[4]]
class_names = [class_names[1],class_names[2],class_names[4]]
rbd_chains = [rbd_chains[1], rbd_chains[2],rbd_chains[4]]
'''

'''
print (ab_names)
print (pdb_names)

pdb_rbd =dict(zip(pdb_names, rbd_chains))

# input

ab_pdb = dict(zip(ab_names, pdb_names))
ab_class = dict(zip(ab_names, class_names))

pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/' for ab in ab_names]))
repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))
fig_dirs = dict(zip(pdb_names,[ '../../output/figs/' + ab +'/' for ab in ab_names]))
merge_dir = '../../output/merge/'
fig_dir = '../../output/figs/'

'''

''' Plot IC50s for antibodies studied '''

def plot_ic50():

    df1 = pd.read_csv('../../data/estimator/NeutSeqData_VH3-53_66.csv')
    df2 = pd.read_csv('../../data/estimator/NeutSeqData_REGN_ngml_units.csv', dtype='object')

    reg = [ab for ab in ab_names_list if 'REGN' in ab]

    df1ss = df1.loc[df1.Name.isin(ab_names_list)]
    df2ss = pd.DataFrame(df2[reg].transpose().iloc[:,0].reset_index().values, columns={'Name','IC50_ngml'})

    df1ss = df1ss[['Name','IC50_ngml']]
    
    df = pd.concat([df1ss, df2ss]).set_index('Name')
    df['IC50_ngml'] = df.IC50_ngml.astype(float)

    sorted = df.sort_values('IC50_ngml')

    plt.subplots(figsize=(1.85,2))
    sb.set(context='paper', style='white') 
    sb.heatmap(data=sorted,cmap='RdBu_r', alpha=0.8)
    plt.xlabel('IC50')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(fig_dir + 'IC50.png', dpi=300)
    plt.show()


''' Plot estimator locations '''

def plot_estimator_locations(ab_name):
    
    estimator = pd.read_csv('../../output/constrain/' + ab_name + '/' + ab_name + '_estimator.csv')
    
    # look at impact of potential mutation at location 

    maxval = estimator.loc[estimator.coefficient > 0].coefficient.max()
    minval = estimator.loc[estimator.coefficient < 0].coefficient.min()

    print (maxval)
    print (minval)
    estpos = estimator.loc[estimator.coefficient > 0].coefficient/maxval 
    estneg = -estimator.loc[estimator.coefficient <= 0].coefficient/minval

    estimator['cf_norm'] = pd.concat([estpos, estneg]).values 
    estimator['pdb_location'] = estimator.pdb_location.astype(int)

    print(estimator.columns)
    print(estimator.pdb_location)

    sorted = estimator.sort_values('pdb_location').loc[abs(estimator.cf_norm) >=1e-3].reset_index(drop=True)

    swt = sorted.loc[sorted.wild_type == True]
    snwt = sorted.loc[sorted.wild_type == False]

    swt['Coefficient'] = swt.cf_norm
    snwt['Coefficient'] = snwt.cf_norm

    pdblocs = swt.pdb_location.append(snwt.pdb_location)
    pdblocs = pdblocs.unique()

    ' FacetGrid is being buggy '
    sb.set(context='paper', style='white')
    plt.figure(figsize=(1.25,2))
    sb.barplot(data=swt,x='mut',y='Coefficient', color='#e31a1c', alpha=0.8)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(fig_dir + 'wt.png', dpi=300)
    plt.show()
    sb.set(context='paper', style='white')
    plt.figure(figsize=(2,2))
    sb.barplot(data=snwt,x='mut',y='Coefficient', color='#1f78b4', alpha=0.8)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(fig_dir + 'notwt.png', dpi=300)
    plt.show()

    return pdblocs 

def plot_energy_locations(ab_name, locs):
    
    energy = pd.read_table('../../output/energy/' + ab_name + '/ddgbind_' + ab_name + '_ab_scanning.txt', sep=',')
    energy['pdb_location'] = energy.pdb_location.astype(int)

    energy['ddg_bind'] = stats.zscore(energy.ddg_bind.values)
    energy = energy.loc[energy['pdb_location'].isin(locs)]
    energy['aa_mut'] = energy.mut.str[-1]

    energy = energy.loc[~energy.aa_mut.isin(['e','o'])]
    print(energy)
    
    colors1 = sb.color_palette("Paired")
    colors2 = sb.color_palette("Spectral")

    colors = colors1 + colors2 
    markers = ['D','o'] 
    plt.figure(figsize=(3.5,7))
    sb.set(context='paper', style='white', font_scale =1.2)
    sb.violinplot(data=energy, y='pdb_location', x='ddg_bind', color='white', scale ='width')
    sb.stripplot(data=energy, y='pdb_location', x='ddg_bind', hue='aa_mut',marker='o',palette=colors, alpha=0.8, s=7,jitter=0.2)
    plt.axhline(0, lw=0.5, color='grey')
    plt.xticks(rotation=90)
    #plt.legend(ncol=7, loc='upper left')
    plt.legend([])
    plt.xlabel('')
    plt.ylabel('$ \Delta \Delta G_{bind}$')
    plt.tight_layout()
    plt.savefig(fig_dir + ab_name + '_ddg.png', dpi=300)
    plt.show()

    
    
    # look at impact of potential mutation at location 

def plot_optimized_ab_locations(ab_name):

    pdb = ab_pdb[ab_name]
    energy_dir = energy_dirs[pdb]


    ' Constrain mutation locations where ddG < 0, wild type == True, based on some cutoffs '
    filename= energy_dir + 'ddgbind_' + ab_name + '_ab_scanning.txt'
    energies = pd.read_csv(filename)

    ' Filter o, a, e, and proline mutations ' 
    energies = energies.loc[~energies.mut_chain.str[-1].isin(['o','a','e','P'])]
    energies['ddg_pos_sign'] = energies.ddg_bind > 0 

    ' Calculate the fitness of each optimization, location most likely to generate better binding '
    print(energies.head())

    grouped = energies.groupby('pdb_location').mean().reset_index()

    sb.set(context='paper', style='white')
    sb.violinplot(data=energies, x='pdb_location', y='ddg_bind', color='white')
    sb.stripplot(data=energies, x='pdb_location', y='ddg_bind', color='red')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
    print (grouped.head())
    exit()

    ' Graph all energies at the locations found by estimator ' 
    cs = energies.sort_values('ddg_bind').reset_index()
    title = 'FoldX binding energies ' + ab_name + ' mutations'
    labels=['',"$\Delta\Delta G_{bind}$"]

    print(cs) 
    exit()

    # normalize between between -1 1 
    locations = [27,28,58,95, 96]
    #cs = cs.loc[cs.pdb_location.isin(locations)]
    cs['norm'] = stats.zscore(cs.ddg_bind.values)

    title = 'FoldX binding energies ' + ab_name + ' all locations'
    pu.violinplot(cs, x='pdb_location', y='norm', hue = 'mut', palette = 'Set1', title= title, figsize=(6,5), labels=labels)

    ap.constrain (constraintype ='energy', constrainfile=filename, antibody=ab_name, cutoff = [-1e4, 0.4], top=1000, out_dir=constrain_dir)
    constrained_energies = pd.read_csv(constrain_dir + ab_name +'_energies.csv')

    print(constrained_energies)
    




''' Calculate the ddg bind of WT antibody and mutations of RBD '''
def viral_scanning_binding():

    for ab in ab_names: 
        print (ab)
        pdb = ab_pdb[ab]
        scan_dir = scan_dirs[pdb]
        energy_dir = energy_dirs[pdb]
        pdb_file = pdb + '_Repair'
        pdb_less_file = pdb + '_Repair_less_ab_Repair'

        ' Find ddG (antibody/receptor) binding after viral scanning '
        ap.energy (antibody=ab, pdb=pdb_file, pdb_less=pdb_less_file, scantype = 'virus', energy_type='ddgbind', indir = scan_dir, outdir=energy_dir)

''' UMAP antibody distance '''
def learn_antibody_distance(antibody_virus_landscape, alg):

    mm = antibody_virus_landscape
    mm = mm.dropna(axis=0)

    ''' Graph distances between antibodies  '''
    # abs = mm.columns.values[1:]
    fitness_data = mm.iloc[:,:].values.transpose()

    if alg == 'UMAP':
        reducer = umap.UMAP(n_neighbors=2,min_dist=1.0,n_components=2)
    if alg == 'MDS':
        reducer = MDS(metric=False, n_components=2)

    #scaled_fitness_data = StandardScaler().fit_transform(fitness_data)
    embedding = pd.DataFrame(reducer.fit_transform(fitness_data), columns={'dim1','dim2'})

    return embedding 




# FITNESS LANDSCAPE VIRUS AND ANTIBODY 
def combine_data():

    merged = pd.DataFrame()
    mergedsens = pd.DataFrame() 

    for ab in ab_names:
        title = 'RBD scanning ' + ab + ' ' +  ab_class[ab]
        pdb = ab_pdb[ab]
        energy_dir = energy_dirs[pdb]

        ddg_label = '$\Delta \Delta G_{bind}($' + ab + '/' + '$\Delta$RBD)' 
        sens_label = '$S$'
        filen ='ddgbind_' + ab + '_virus_scanning.txt'
        df = pd.read_table(energy_dir + filen, ',')
        df['antibody'] = ab
        df['mut_nochain'] = df.mut.str[0] + df.mut.str[2:]

        # Filter o, a, e, and proline mutations from FoldX
        df['mut'] = df.mut.str[-1]
        
        lmut = ['o','a','e','P']
        df = df[~df.mut.isin(lmut)]

        # consider only positions 400 to 520 
        df = df.loc[(df.pdb_location >=400) & (df.pdb_location <=520)]

        # Normalize the ddg_bind
        minval = df.ddg_bind.min()
        maxval = df.ddg_bind.max()

        dfn = df.loc[df.ddg_bind <=0]
        dfn['ddg_bind_norm'] = -dfn['ddg_bind']/minval

        dfp= df.loc[df.ddg_bind >0]
        dfp['ddg_bind_norm'] = dfp['ddg_bind']/maxval
        df = pd.concat([dfn, dfp])

        #pdb_file_name = pdb + '.pdb'
        #labeled_chains = structure.label_chains(pdb)
        #ep = structure.find_epitopes(pdb_dirs[pdb], pdb_file_name, labeled_chains, 4)    

        #epitope = pd.read_csv('../../data/location/' + pdb + '_epitopes.csv')        
        #proximity = epitope.number_virus.unique()

        # remove chain and mutation that end with e or o

        locations = pd.DataFrame(df.pdb_location.unique(), columns={'pdb_location'})

        print(df.head())
        sensitivity = df.loc[df.ddg_bind_norm >0]

        sensitivity = sensitivity.groupby(by=['antibody', 'pdb_location']).count().reset_index()

        print(sensitivity.head())
        print(locations.head())
        msens = sensitivity.merge(locations, on = 'pdb_location', how='right')
        msens = msens.fillna(0)
        msens['antibody'] = ab
        msens['S'] = msens.ddg_bind_norm
        sensitivity = msens[['pdb_location','S']].set_index('pdb_location')

        mergedsens = pd.concat([mergedsens, msens[['pdb_location', 'S', 'antibody']]])
        #vlines = locations.loc[locations['pdb_location'].isin(proximity)].index

        if plot: 
            sb.set(context='paper', style='ticks')
            f, axes = plt.subplots(1, 1, figsize=(10, 2), sharex=True, sharey= False)
            #ax0 = sb.heatmap(data=sensitivity.transpose(), vmin = 0, vmax = 21, ax=axes[0], square=True, cmap='RdBu_r', cbar = False, alpha = 0.8)
            #ax0.set_ylabel('')
            #ax0.set_xlabel('')
            sb.violinplot(data=df, x='pdb_location', y='ddg_bind_norm',color='white',scale='width', ax=axes,inner = None, lw=0.5)
            ax2= sb.stripplot(data=df, x='pdb_location', y='ddg_bind_norm', color='red', alpha=0.3, s=2, ax=axes)
            #plt.vlines(vlines, ymin = -1.5, ymax= 1.5, color='#ffeda0', alpha=0.6, lw=2)
            plt.xticks(rotation=90)
            plt.xlabel('RBD location')
            plt.subplots_adjust()
            f.suptitle(title)
            plt.xlim([0,120])
            plt.ylim([-1.5,1.5])
            plt.xlabel('')
            plt.ylabel(ddg_label)
            plt.tight_layout()
            plt.savefig(fig_dir + title + '.png', dpi =dpi)
            plt.show()
        
        merged = pd.concat([df, merged])

    # group antibodies then calculate virus sensitivities  # of positive ddgs 
    mergedsens.to_csv(merge_dir + 'ddgbind_virus_scanning_sensitivity.csv')

    merged = merged[['antibody',  'mut_nochain', 'ddg_bind_norm']]
    merged = merged.dropna()

    merged.to_csv(merge_dir + 'ddgbind_virus_scanning.csv')

    mergedsens['S'] = mergedsens['S']/20

    meltedmergedsens = mergedsens.pivot(index= ['pdb_location'] ,columns=[ 'antibody'], values =['S'])
    meltedmergedsens = meltedmergedsens.fillna(0)

    meltedmerged = merged.pivot_table(index='mut_nochain',columns='antibody')['ddg_bind_norm']
    meltedmerged = meltedmerged.fillna(0)

    print(meltedmerged.head())

    # Filter by location 
    meltedmergedsens = meltedmergedsens.loc[meltedmergedsens.index>=400]
    meltedmergedsens = meltedmergedsens.loc[meltedmergedsens.index < 520]

    # Filter by non WT 
    meltedmerged['WT'] = meltedmerged.index.str[0] == meltedmerged.index.str[-1]

    meltedmerged = meltedmerged.loc[meltedmerged.WT== False]
    meltedmerged = meltedmerged.drop(['WT'], axis=1)

    # Sort by location 
    meltedmerged['vloc'] = meltedmerged.index.str[1:-1].astype(int)

    sorted = meltedmerged.sort_values('vloc')
    vloc = sorted.vloc.values
    sorted = sorted.drop('vloc',axis=1)

    '''
    plt.figure(figsize=(10,2.5))
    sb.heatmap(data=sorted.transpose(),square=False, cmap='RdBu_r',xticklabels= 40,alpha=2.5, cbar=True)#, col_cluster=False)
    plt.tight_layout()
    plt.savefig(fig_dir + 'virus_scan_WT_ab_hm.png', dpi =dpi)
    plt.show()
'''
    # add the RBD importances here 
    tol = 0.003
    rbd_imp = pd.read_csv('../../output/estimator/single_mut_effects_cleaned_coefficients.csv')

    rbd = rbd_imp.loc[rbd_imp.importances >= tol]
    rbd['location'] = rbd['Unnamed: 0'].str[1:].astype(int)

    rbd.to_csv(merge_dir + 'tmp.csv')
    m = pd.merge(rbd,meltedmergedsens, how='right', left_on='location', right_on='pdb_location')

    # run MDS analysis or UMAP analysis  
    mm = sorted.reset_index()

    print (mm.head())
    print (mm.columns)
    mm['gloc'] = mm.mut_nochain.str[1:-1].astype(int)

    # meman, mean of postives, 
    grouped = mm.groupby('gloc').mean().reset_index().iloc[:,1:]
    grouped['gloc'] = mm.gloc.unique()

    pos = grouped.iloc[0:-1] > 0 

    print(pos.head())

    sb.set(context='paper', style='ticks')
    f, axes = plt.subplots(1, 1, figsize=(10, 2), sharex=True, sharey= False)
    #sb.barplot(data=grouped, x= 'gloc', y='B38', alpha=0.4)
    sb.lineplot(data=grouped, x= 'gloc', y='B38')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

    print(grouped.head())
    cols = mm.columns[1:-1]
    print(cols)
    
    # plot antibody distance
    umap_antibody_distance(mm)

    # plot RBD 


#plot_ic50()        

#ab_name = 'C105'

#pdblocs = plot_estimator_locations(ab_name)


#pdblocs = [27, 28, 58, 95, 96,]
#plot_energy_locations(ab_name, pdblocs)
#plot_optimized_ab_locations(ab_name)


#viral_scanning_binding()
#exit()
#combine_data()
