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

ab_names = ['CB6', 'C105', 'B38','CC12.1', 'ACE2']
pdb_names =['7c01','6XCM', '7bz5','6xc3', '6m0j']
has_mutations =[True, True, False, False]

ab_muts = dict(zip(ab_names, has_mutations))
ab_pdb = dict(zip(ab_names, pdb_names))




def get_subset_ddg(melted, yval = 'ddg', cutoff= 20):

    melted = melted.loc[melted[yval]>=0]
    melted.to_csv(out_tmp + 'tmp.csv')
    y= melted.sort_values(by='ddg', ascending=False)
    return y 

def get_cutoff_ddg(melted, yval = 'ddg'):
    cutoff = 3 
    m = melted.loc[melted[yval]>=cutoff]
    y= m.sort_values(by='ddg', ascending=False)
    return y


'''
# mutations of RBD in context of C105_wt

# mutations of RBD in context of C105_mt  
p1 = '6XCM_Repair_TH28D_Repair'
p2 = '6XCM_Repair_TH28D_Repair_less_ab_Repair' #less ab, mutate virus 

#plot_ddg_distribution(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')
#plot_ddg_heatmap(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')

ab_name = 'C105'
pdb_name = '6XCM_Repair'
less_str = '_less_ab_Repair'
'''

def get_all_ddgs_for_antibody_mutations(ab_name): 
    mut_filename = ab_name + '_ab_mutations.txt'
    muts = pd.read_table(ab.pdb_path + ab_name + '/' + mut_filename, header=None)
    muts = muts.iloc[:,0].str.rstrip(';')
    
    ddg_array = []
    for mut in muts:
        mut_pdb = pdb_name + '_' + mut + '_Repair'
        mut_name = pdb_name + '_' + mut
        ddg = ab.calculate_ddg_bind(ab_name, mut_pdb , mut_pdb + less_str, mut_name = mut_name) 
        ddg_array.append(ddg)
    return ddg_array

'''
# WT C105 antibody
ab_name = 'C105'
pdb_name = '6XCM_Repair'
less_str = '_less_ab_Repair'

mut_pdb = pdb_name
mut_name = pdb_name
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
ddgs = ddgs.drop(['sumval'], axis=1)

'''

def heatmap_ab_sensitivity(ddgs,title='hm'):

    ddgs = ddgs.transpose()

    sb.set(context=context, style=style)
    plt.figure(figsize=(12,2.5))
    sb.heatmap(ddgs, square=False,cmap='RdBu_r')
    plt.title('ddg(ab/RBD)<0 normalized')
    plt.tight_layout()
    plt.savefig(out_fig + title + '.png', dpi=dpi)
    plt.show()


# scatterplot facet grid of this data wt vs others 

'''
wt = 'ddg_6XCM_Repair'

muts_C105 = ddgs.columns[1:6]
print (muts_C105)

# subtract wt from mutant: MUT-WT this will give positions where it is less or more binding 
#ddgs = ddgs - ddgs[wt]

dgwtmut = pd.DataFrame()

dgwtmut[wt+'LESS'] = ddgs[wt]- ddgs[wt]
for m in muts_C105:
    col = m +'_LESS' 
    dgwtmut[col] = ddgs[m]- ddgs[wt]

#heatmap_ab_sensitivity(ddgs,title='all_abs__')
#heatmap_ab_sensitivity(dgwtmut,title='C105_less_hm_')

'''



# normalize data between min and max 
def normalize_data(data, colname='coefficient', min = '-1', max ='1'):

    dss = data
    dssn = dss.loc[dss[colname] <= 0]
    dssp = dss.loc[dss[colname] > 0]
    dssn['norm_'+colname] = -dssn[colname]/dssn[colname].min()
    dssp['norm_' +colname] = dssp[colname]/dssp[colname].max()
    y = pd.concat([dssn, dssp])
    yval = 'norm_'+colname
    y = y.sort_values(by=yval, ascending=False)

    return y, yval 


def get_estimator (estimator_file, ab_name,  filtername, filterval, cutoffmin=0, cutoffmax=0):

    # cutoff 1e-10 for nussenzweig and 1e-5 for vh3-53/66                                              
    estimator = pd.read_csv(est_path + estimator_file)
    estimator['fasta_location'] = estimator.location.astype(str)
    estimator = estimator.loc[(estimator.Name == ab_name) & (estimator[filtername] == filterval)]
    estimator = estimator.loc[(estimator.coefficient < cutoffmin ) | (estimator.coefficient > cutoffmax)]

    return estimator 

def compare_estimator_ddgs(estimator_file,ab_name, pdb_name, scantype='virus'):

    estimator = get_estimator (estimator_file, ab_name,  filtername='chain', filterval='H', cutoffmin=-1e5, cutoffmax=1e5)
    coeff, coeffval = normalize_data(estimator,colname='coefficient')    
    plot_scatter(ab_name,coeff, coeffval) 

    less_virus = '_less_virus_Repair'
    ddgs = ab.calculate_ddg_bind(ab_name,pdb_name, pdb_name + less_virus, scantype='ab', mut_name='')
    ddgs['location'] = ddgs.mut.str[3:-1]
    ddgs['mut'] = ddgs.mut.str[-1:]
    ddgs, ddgsval = normalize_data(ddgs, colname='ddg_')
    plot_scatter(ab_name, ddgs, ddgsval) 

    plt.figure(figsize=(12,4))
    plot_ddg_stripplot(ab_name, ddgs,x='location', y=ddgsval, hue='mut', title=ab_name + 'ab_scann_ddg', subset=False)

    coeff['loc_mut'] = coeff.location.astype(str)  + coeff.aa.astype(str)
    ddgs['loc_mut'] = ddgs.location.astype(str)  + ddgs.mut.astype(str)

    merged  = ddgs.merge(coeff, on='loc_mut')

    merged['pdb_location'] = merged.location_x.values.astype(object)
    merged = merged[['loc_mut', ddgsval, coeffval, 'location_x', 'chain']]
    merged['ab_name'] = ab_name
    merged.to_csv(out_tmp + 'tmp.csv')

    plot_scatter_joint(ab_name, data=merged, x=coeffval, y=ddgsval, hue='loc_mut') 
    
    return merged

def compare_all_estimator_ddgs():

    est_file =  'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv' 
    m1 = compare_estimator_ddgs(est_file,ab_name='B38', pdb_name ='7bz5_Repair', scantype='ab')
    # compare_estimator_ddgs(ab_name='CC12.1', pdb_name ='6xc3_Repair', scantype='ab')
    m3 = compare_estimator_ddgs(est_file, ab_name='CB6', pdb_name ='7c01_Repair', scantype='ab')
    m4 = compare_estimator_ddgs(est_file, ab_name='C105', pdb_name ='6XCM_Repair', scantype='ab')
    m5 = compare_estimator_ddgs(est_file, ab_name='CV30', pdb_name ='6xe1_Repair', scantype='ab')
    
    m = pd.concat([m1,m3, m4, m5])
    
    ddgsval = 'norm_ddg_'
    coeffval = 'norm_coefficient'
    plot_scatter_joint(ab_name='all', data=m, x=coeffval, y=ddgsval, hue='ab_name', leg=True)

def convert_to_simple_mutation(ddgs):

    ddgs['wt_three'] = ddgs.mut_chain.str[0:3]
    ddgs['location'] = ddgs.mut_chain.str[4:-1]
    ddgs['mut'] = ddgs.mut_chain.str[-1:]

    u.write(ddgs)
    ddgs = ddgs.loc[ddgs.wt_three !='H2S']


    wtone = [ab.aa_3to1_dic[wt] for wt in ddgs.wt_three.values]

    ddgs['wt'] = wtone
    ddgs['mutation'] = ddgs.wt + ddgs.location + ddgs.mut
    return ddgs 

def compare_bloom_ddgs():
    f = 'single_mut_effects_cleaned.csv'

    bloom = pd.read_csv(ml_path + f)
    df = bloom [['mutation', 'bind_avg']]

    bloom, bloomval = normalize_data(df,colname='bind_avg')    

    plot_scatter('bloom', bloom, bloomval)

    ab_name = 'ACE2'
    pdb_name = '6m0j_Repair'
    less_ab = '_less_ACE2_Repair'

    df = ab.calculate_ddg_bind(ab_name,pdb_name, pdb_name + less_ab, scantype='virus', mut_name='')
    ddgs = convert_to_simple_mutation(df)
    ddgs, ddgsval = normalize_data(df,colname='ddg_')    

    plot_scatter(ab_name, ddgs, ddgsval)

    merged = bloom.merge(ddgs, on='mutation') 

    merged.to_csv(out_tmp + 'tmp.csv')
    bloomval = 'norm_bind_avg'
    ddgsval = 'ddg_'

    plot_scatter_joint(ab_name, merged, x=bloomval, y=ddgsval,hue='mutation', leg=False)


 
def get_bloom_data(normalize=True):
    f = 'single_mut_effects_cleaned.csv'
    bloom = pd.read_csv(ml_path + f)
    df = bloom [['mutation', 'bind_avg']]
    df['bind_avg_bloom'] = -df.bind_avg.values
    bloom = df 
    bloomval = 'bind_avg_bloom'
    if normalize: 
        bloom, bloomval = normalize_data(df,colname='bind_avg_bloom')
    return bloom, bloomval


def get_ab_ddg_data(abname, mutation = None, scantype='virus', normalize = True):

    dgpath = data_path + 'ddg/' + abname + '/'
    pdbname = ab_pdb[abname]

    mutstr , ddglabel, ablabel  = '', '', abname
    if mutation != None: 
        mutstr = '_' + mut + '_Repair'
        ab_label = abname + '_' + mut

    f1 = 'PS_' + pdbname + '_Repair' + mutstr + '_' + scantype + '_scanning_output.txt'
    f2 = 'PS_' + pdbname + '_Repair' + mutstr + '_less_ab_Repair_' + scantype + '_scanning_output.txt'
    t1 = pd.read_table(dgpath + f1)
    t2 = pd.read_table(dgpath + f2)
    ddg = ab.calculate_ddg_bind(abname,pdbname +'_Repair', pdbname +'_Repair_less_ab_Repair', scantype=scantype)
    ddglabel = 'ddg'

    if normalize: 
        ddg, ddglabel = normalize_data(ddg, colname='ddg_')

    ddg['abname'] = ablabel 
    ddg = convert_to_simple_mutation(ddg)
    return ddg



def create_ab_fitness_matrix():
    ddgall = pd.DataFrame()

    for ab_name in ab_names:
        path_dg = data_path + 'ddg/' + ab_name + '/'
        path_mut = data_path + 'pdb/' + ab_name + '/'

        merged, bloomval, ddgval = merge_ab_bloom(ab_name, mut ='', normalize=True)
        ddgall = pd.concat([merged, ddgall])

        if ab_muts[ab_name] == True:
            print(ab_name)
            mutations = pd.read_table(path_mut + ab_pdb[ab_name] + '_ab_mutations.txt', header=None)
            for mut in mutations.iloc[:,0]:
                print (mut)
                merged, bloomval, ddgval = merge_ab_bloom(ab_name, mut=mut.strip(';'), normalize=True)
                ddgall = pd.concat([merged, ddgall])
    
                
    ddgall['fit'] = ddgall[bloomval] - ddgall[ddgval]
    suff = ''
    ablist = ['C105_', 'CB6_', 'CC12.1_', 'B38_']

    ddgall= ddgall.loc[ddgall.ab_name.isin(ablist)]

    ddgall = ddgall.pivot(index='mutation', columns='ab_name')['fit']    
    ddgall = ddgall.dropna(axis=0)
    ddgall.to_csv(out_tmp + 'matrix.csv')

    sb.set(context='paper')
    sb.clustermap(ddgall.transpose(),cmap='RdBu', figsize=(12,6))
    plt.tight_layout()
    plt.savefig(out_fig +'hm_' + suff + 'fitness.png', dpi = dpi)
    plt.show()

    sign_mat = pd.DataFrame(np.sign(ddgall.values), columns=ddgall.columns, index=ddgall.index)
    sign_mat.to_csv(out_tmp + 'sign_matrix.csv')

    tmat = (sign_mat == -1 )

    # if all are False then the problem will be infeasible
    tmat['sum_false'] = tmat.sum(axis=1)
    not_covered = tmat.loc[tmat.sum_false == 22]
    not_covered  = not_covered.reset_index()
    not_covered['location']  = not_covered.mutation.str[1:-1]

'''

# S2: get antibody estimators and ddgs 
estimator_file =  'NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv' 
ab_name = 'CB6'
pdb = ab_pdb[ab_name] + '_Repair'
pdb_less = pdb + '_less_virus_Repair'
estimator = get_estimator (estimator_file, ab_name,  filtername='chain', filterval='H', cutoffmin=-1e-5, cutoffmax=1e-5)
ddgs = ab.calculate_ddg_bind(ab_name,pdb, pdb_less, scantype='ab', mut_name='')

# add to structure 
mutable_estimator = ab.get_antibody_mutation_positions(estimator, cutoffmin =-1e-5, cutoffmax=1e-5, topmuts=10)
u.write(mutable_estimator)

mutable_locations = mutable_estimator.pdb_location.unique()
ddg_estimator = ddgs.loc[ddgs.pdb_location.isin(mutable_locations)]
mutation_e = ddg_estimator.mut.str[-1] =='e' 
mutation_o = ddg_estimator.mut.str[-1] =='o' 
muts = mutation_e | mutation_o
ddg_estimator['muts'] = muts.values

u.write(ddg_estimator)
ddg_estimator = ddg_estimator.loc[ddg_estimator.muts == False]
ddg_estimator['is_mutable'] = ddg_estimator['ddg_bind'] < 0.8

topmuts = 3
all_mutable = pd.DataFrame()

for wt in ddg_estimator.wt.values: 

    ddgss = ddg_estimator.loc[ddg_estimator.wt == wt]
    top = ddgss.loc[(ddgss.ddg_bind < 0.8) & (ddgss.wildtype == False)]
    topthree = top.loc[top.ddg_bind.astype(float).nsmallest(n=topmuts).index]

    all_mutable = pd.concat([all_mutable, topthree])
    ddgss.sort_values(by='ddg_bind', inplace=True)
    plt.figure(figsize=(3,2.5))
    sb.set(context='paper', style='ticks')
    sb.barplot(data= ddgss, x='mut', y="ddg_bind",hue = 'is_mutable', palette = 'RdBu', dodge=False)
    plt.gca().legend().set_title('')
    plt.title(ab_name + ' ' + wt)
    plt.axhline(y=0.8, ls='--', lw=0.5, color='gray')
    plt.axhline(y=0.0, ls='-', lw=0.5)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_fig + ab_name +'_' +  wt + '.png')
    plt.close()
    
all_mutable = all_mutable.drop_duplicates()
all_mutable['antibody'] = ab_name
all_mutable['pdb'] = ab_pdb[ab_name]
sorted = all_mutable.sort_values('ddg_bind')
u.write(all_mutable)

plt.figure(figsize=(3,2.5))
sb.set(context='paper', style='ticks')
sb.barplot(data= sorted, x='mut', y="ddg_bind", color = '#3182bd', dodge=False, alpha=0.5)
#plt.gca().legend().set_title('')
plt.axhline(y=0.0, ls='-', lw=0.5)
plt.xticks(rotation=90)
plt.title(ab_name + ' mutations') 
plt.tight_layout()
plt.savefig(out_fig + ab_name + '_all_mutable.png')
plt.close()
exit()


sorted = ddg_estimator.sort_values(by ='ddg_bind')
bar_plot_ddgs(data=sorted, ab_name=ab_name, figsize = (13,3))

# get top 10 and  lower 10 
topmuts = 15 
top = sorted.loc[sorted.ddg_bind.astype(float).nlargest(n=topmuts).index]
bottom = sorted.loc[sorted.ddg_bind.astype(float).nsmallest(n=topmuts).index]

# top = sorted.loc[sorted.ddg_bind.astype(float) >= 1]
# bottom = sorted.loc[sorted.ddg_bind.astype(float) <= -0.5]

sorted_mutations = pd.concat([top, bottom]).sort_values(by ='ddg_bind')
bar_plot_ddgs(data=sorted_mutations, ab_name=ab_name, figsize = (6,3))
'''
