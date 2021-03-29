import numpy as np
import cvxpy as cp
import pandas as pd
import seaborn as sb
import argparse
import matplotlib.pyplot as plt
import os 
import utils 
import colors as fmt 
import manuscript_analysis as ma 
from itertools import combinations
import scipy.stats as stats

def setup_parser():
    parser = argparse.ArgumentParser(prog = 'OPTIMIZE', description = 'Uses convex combinatorial optimization to compute optimal antibody cocktail for provided virus mutants.')
    parser.add_argument('--filepath', type = str, required = True, dest = 'filepath', help = 'Filepath of fitness matrix dataframe', metavar = 'FILEP')
    parser.add_argument('--filename', type = str, required = True, dest = 'filename', help = 'Filename of fitness matrix dataframe (do not include file type extension)', metavar = 'FILEN')
    parser.add_argument('--virus_filepath', type = str, required = True, dest = 'v_filepath', help = 'Filepath of fitness matrix dataframe for virus mutants', metavar = 'VFILEP')
    parser.add_argument('--virus_filename', type = str, required = True, dest = 'v_filename', help = 'Filename of fitness matrix dataframe (do not include file type extension) for virus mutants', metavar = 'VFILEN')
    parser.add_argument('--p', type = float, required = True, dest = 'p_unmanaged', help = 'Maximum proportion of viruses not covered by the antibody cocktail', metavar = 'P_UNMANAGED')
    parser.add_argument('--g1', type = float, required = True, dest = 'gamma1', help = 'Gamma 1 value to use for penalization of number of antibodies chosen', metavar = 'GAM1')
    parser.add_argument('--g2', type = float, required = True, dest = 'gamma2', help = 'Gamma 2 value to use for weighting infectivity of virus mutants', metavar = 'GAM2')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser


def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/cocktail/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')



''' Private functions '''
def import_fitness_landscapes(antibody_fitness_file, virus_fitness_file):
    
    ''' Import fitness landscapes:  ligand1 -> mutations receptor1,  ligand2 -> mutations receptor 1  
        eg: antibody-> mutations RBD,  ACE2 -> mutations RBD 
        Returns dataframes ligand 1 and ligand 2 fitness landscapes 
    '''

    cocktail_path = '../../output/cocktail/'
    
    abdf = pd.read_csv(antibody_fitness_file, sep=',', header=0, index_col=0)

    ''' Drop C135 because of missing structure ''' 
    abdf = abdf.drop('C135', axis=1)
    vdf = pd.read_csv(virus_fitness_file, sep=',', header=0, index_col=0)

    merged = abdf.merge(vdf,how='inner', on='mut')
    merged = merged.dropna(axis=0)
    virus_fitness = merged[['ACE2']]
    antibody_fitness = merged.drop('ACE2',axis=1)    

    return antibody_fitness, virus_fitness


def compute_cocktail(antibody_fitness,virus_fitness, k, lmbda1, lmbda2, noise):

    ''' Compute antibody cocktail given k, l1, l2 '''

    fitness_matrix = antibody_fitness.values + noise

    virus_matrix = virus_fitness.values
    m,n = fitness_matrix.shape
    num_unmanaged = int(m*k)
    c = cp.Variable(n, boolean = True)

    incidence_matrix = np.sign(-fitness_matrix).clip(min=0) # less than zero 
    virus_incidence_matrix = np.sign(virus_matrix -fmt.MIN_ACE2_DDG_BIND).clip(min=0) # ones that greater than zero clip to zero
    constraints = [cp.sum(c) >= 1, cp.sum_smallest((incidence_matrix@c), num_unmanaged+1) >= 1]

    # scale between [-1,1] for each antibody, here we measure the robustness of each Ab wrt mutations of virus  
    '''
    print(antibody_fitness.head())
    print(antibody_fitness.max())
    print(antibody_fitness.min())
    print((antibody_fitness>0).sum())
    print((antibody_fitness.C105_TH28I_YH58F ==  antibody_fitness.C105).sum())
    '''

    # objective = cp.Minimize(lmbda1*cp.norm1(c)+cp.matmul(cp.sum(fitness_matrix, axis=0), c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    # mean of antibody cocktail: 
    #objective = cp.Minimize(lmbda1*cp.norm1(c)+cp.matmul(cp.sum(fitness_matrix, axis=0), c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    #objective = cp.Minimize(lmbda1*cp.norm1(c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    objective = cp.Minimize(lmbda1*cp.norm1(c)+ cp.matmul(cp.sum(fitness_matrix, axis=0),c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value


''' API '''
def run_simulations(antibody_fitness, virus_fitness, coverage,gamma1, gamma2, noise_sims=0):

    ''' Runobjective = cp.Minimize(lmbda1*cp.norm1(c)+cp.matmul(cp.sum(fitness_matrix, axis=0), c)-lmbda2*incidence_matrix@c@virus_incidence_matrix)
    problem = cp.Problem(objective, constraints)
    result = problem.solve() simulations for different mutation coverages and regularization parameters 
        antibody_fitness: file with fitness landscape of antibodies with respect to virus binding
        virus_fitness: file with fitness landscape of virus fitness with respect to cell receptor binding  

        coverage: float or array of minimum virus coverage
        gamma1: float or array  of lambda1 values, tuning parameter controlling the  number of antibodies in mix 
        gamma2: float, or array of lambda 2 values, tuning parameter controlling the infection  
    '''    

    abdf, vdf = import_fitness_landscapes(antibody_fitness, virus_fitness)

    allresults = pd.DataFrame()
    allresults['antibody'] = abdf.columns

    ks = coverage 

    allsims = pd.DataFrame()

    ''' Subject this to noise ''' 

    for i in range(noise_sims+1): 

        noise_str = 'noise_'+ str(i)

        ''' Generate random noise ''' 
        random_noise = np.zeros(abdf.shape)
        if i > 0:
            random_noise = np.random.normal(loc=0.0, scale=1, size=abdf.shape)*1e-2

        ''' Loop through coverage '''
        allcov = pd.DataFrame()
        for k in ks:
            print('coverage=' + str(k) + noise_str) 

            '''Grid search on gamma1, gamma2 '''
            for g1 in gamma1: 
                for g2 in gamma2:
                    results = compute_cocktail(abdf,vdf, k, g1, g2, random_noise)
                    choices = abdf.columns[results.astype(bool)]
                    colname = str(1-k) + '_' + str(g1) + '_' + str(g2)
                    allresults[colname]=results

            ''' Find minimum abs ''' 
            abt = allresults.set_index('antibody')
            abt = abt.transpose() 
            abt['num_abs']  = abt.sum(axis=1)
            tmp= pd.Series(abt.index).str.split(pat='_',  expand=True)

            abt ['cov'] = tmp.iloc[:,0].values
            abt ['gamma1'] = tmp.iloc[:,1].values
            abt ['gamma2'] = tmp.iloc[:,2].values
            abt ['noise'] = noise_str

            abt.to_csv('../../output/cocktail/cocktails_'+ str(k) + noise_str + '.csv', index=False)

        allsims = pd.concat([allsims, abt]) 

    allsims.to_csv('../../output/cocktail/cocktails_allsims.csv', index=False)
        
    return abt

def graph_simulation(cocktail_file):

    allres = pd.read_csv(cocktail_file)

    ab_class_dict, ab_gene_dict,ab_label_dict = utils.get_antibody_properties()
    ablabels = [ab_label_dict[ab] for ab in allres.antibody.values]
    replace_dict = dict(zip(allres.antibody.values,ablabels))
    allres = allres.replace(replace_dict)
    allres= allres.set_index('antibody')

    figsize = (4.5,4.5)    
    plt.figure(figsize=figsize)
    sb.set(context='paper', style='white')
    sb.heatmap(data=allres ,cmap=colors.INCOCKTAIL,square=True, linecolor='white', linewidths=0.5, cbar=None, alpha=0.7)
    plt.xticks(rotation = 90)
    plt.xlabel('coverage_g1_g2')
    plt.tight_layout() 
    plt.savefig('../../output/manuscript/combinations_' +cocktail_file[-7:-4] + '.png', dpi=300)


def calculate_cocktail_distance():

    figname = 'graph_cocktail_distance'

    data = pd.read_csv('../../output/cocktail/cocktail_virus_fitness.csv')

    data['location'] = data.mut.str[1:-1]
    data['uncovered'] = ~data.iscovered.values 
    data['sortcol'] = data.cocktail.str[1:]
    
    ''' Merge with ACE2 binding ''' 
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)[['bind_avg','mutation']]
    merged = data.merge(ace2df, how='inner', left_on='mut', right_on='mutation')
    merged = merged.drop('mutation', axis=1)
    merged = merged.loc[merged.numabs>2] # any other combinations have less than two cocktails to compare

    allstats= pd.DataFrame()
    for numab in merged.numabs.unique():
        print('Number antibodies '+ str(numab))
        filtered = merged.loc[merged.numabs == numab]
        print(filtered)
        cocktails = filtered.cocktail.unique()

        print(cocktails)
        combs = list(combinations(cocktails, 2))
        print(combs)
        ''' Perform t test and compare accross locations designed ab vs wt  '''
        test= 'ttest'
        sig_min=5e-2
        fc_min=0
        
        for comb in combs: 
            stat_test = pd.DataFrame()
            for loc in filtered['location'].unique():
                filtered_loc = filtered.loc[filtered['location'] == loc]
                filtered_loc = filtered_loc[['mut','location','cocktail','bind_avg','max_neutralization']]                    
                
                c1= filtered_loc.loc[filtered_loc.cocktail  == comb[0]].max_neutralization
                c2= filtered_loc.loc[filtered_loc.cocktail  == comb[1]].max_neutralization

                c1mean = c1.mean()
                c2mean = c2.mean()
                    
                res = stats.ttest_ind(c1, c2)          
                stat_test = stat_test.append(pd.Series([loc, res.pvalue, c1mean, c2mean]), ignore_index=True)

        
            ''' Extract statistically significant locations '''                    
            stat_test = stat_test.rename(columns={0:'location', 1:'pval', 2:comb[0], 3:comb[1]})
            stat_test['significant'] = stat_test.pval < sig_min
            stat_test['fold_change'] = np.abs(stat_test[comb[0]]/stat_test[comb[1]])

            melted = stat_test.melt(value_vars=[comb[0],comb[1]], value_name='mean', id_vars=['location','significant','pval'])
            melted['numab'] = numab
            allstats = pd.concat([melted,allstats], axis=0)

    
        allstats = allstats.loc[allstats['significant'] ==True]
        allstats.to_csv('../../output/cocktail/cocktail_differences.csv', index=False)


def graph_robustness():

    figname = 'graph_cocktail_distance'

    ab_class_dict, ab_gene_dict,ab_label_dict = utils.get_antibody_properties()

    allstats = pd.read_csv('../../output/cocktail/cocktail_differences.csv')
    cv = pd.read_csv('../../output/cocktail/cocktail_virus_fitness.csv')
    neut = pd.read_csv('../../output/cocktail/cocktail_neut_fitness.csv')

    '''ACE2 binding'''
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)[['bind_avg','mutation']]
    abfit = utils.get_antibody_fitness_landscape().merge(ace2df, how='inner', left_on='mut', right_on='mutation')
 
    cocktails = ['C' + str(i) for i in range(3,18)]
    cocktails = cocktails[0:1] + cocktails[5:7] + cocktails[10:11]+ cocktails[12:15]

    abcolors = fmt.CLASS_AB_COLOR_DICT

    ''' Max neutralization combined with ACE2 binding ''' 
    cvace = cv.merge(ace2df, how='inner', left_on='mut', right_on='mutation')
    cond = (cvace.max_neutralization > fmt.MAX_AB_DDG_BIND) & (cvace.bind_avg >= fmt.MIN_ACE2_DDG_BIND)
    escape = cvace.loc[cond]

    grouped = escape.groupby(['cocktail']).count().reset_index()
    grouped  = grouped.loc[grouped.cocktail.isin(cocktails)]
    grouped = grouped[['cocktail', 'mut']]
    grouped['noise'] = 0 

    escapedf = pd.DataFrame()
    escapemuts = pd.DataFrame()

    numsims = 10
    escapevars = []

    ''' Generate random noise ''' 
    random_noise = np.random.normal(loc=0.0, scale=1, size=[numsims,len(abfit.mut)])*1e-2

    sall = pd.DataFrame()
    for j,c in enumerate(cocktails):

        ss = cvace.loc[cvace.cocktail == c][['mut','cocktail','max_neutralization','bind_avg']]

        colnames = ['p' + str(i+1) for i in range(numsims)]
        perturbed = pd.DataFrame(random_noise + ss.max_neutralization.values).transpose()

        perturbed = perturbed.rename(columns=dict(zip(range(numsims), colnames)))
        perturbed['mut'] = ss.mut.values
        ss = ss.merge(perturbed, how='inner')        
        sall = pd.concat([ss, sall])

    cond1 = (sall.bind_avg >= fmt.MIN_ACE2_DDG_BIND)
    cond2 = (sall[sall.columns[4:]] > fmt.MAX_AB_DDG_BIND)        
    cond12 = pd.DataFrame()

    for col  in cond2.columns:
        cond12[col] = cond1 & cond2[col]

    print(cond12.head())
    cond12 = cond12.astype(int)
    cond12['mut']  = sall.mut.values
    cond12['cocktail']  = sall.cocktail.values
    cond12['location']  = sall.mut.str[1:-1]
    print(cond12.head())

    d = pd.DataFrame()
    for c in cocktails: 
        dd = cond12.loc[cond12.cocktail == c]
        print(dd[dd.columns[0:numsims]].sum().values)
        variants = dd[dd.columns[0:numsims]].sum().values
        d[c] = variants

    
    d['simnum'] = range(1,numsims+1)

    print(d)
    melt = d.melt(id_vars='simnum',value_vars=d.columns[:-1]) 



    exit()
    ''' Group by cocktail'''
    g = cond12.groupby('cocktail').count().reset_index()[['location','mut']]
    g['noise'] = j 
    g['cocktail'] = c 
    
    print(g)
    exit()
            
    '''
    s1 = sss['additive-perturbed'].count()
    escape_vars.append(s1)
    
    tmp = pd.DataFrame(escape_vars).transpose()    
    escapedf = pd.concat([escapedf, tmp])

    escapedf = escapedf.rename(columns = dict(zip(range(len(cocktails)),cocktails)))
    melt = escapedf.melt(value_vars=cocktails)

    print(escapemuts)

    for c in cocktails: 
        s = escapemuts.loc[escapemuts.cocktail  == c]
        print (s.head())
        pivot = pd.pivot(data =s , columns ='location', values='mut')
        print(pivot.head())
        exit()

    exit()
    '''
    ''' Now look at locations that are significantly different ''' 


    area  = g.mut.values 

    plt.figure(figsize=(10,2.5))
    sb.scatterplot(data=g, x= 'location', y = 'mut', hue ='cocktail', hue_order=cocktails, palette='RdPu' )
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/scatter_robustness_legend.png', dpi=300, transparent=True)
    plt.show()

    min = escapedf.min(axis=0)
    max = escapedf.max(axis=0)
    vardf = escapedf.var(axis=0)
    mean = escapedf.mean(axis=0)

    plt.figure(figsize=(6,2.5))
    sb.lineplot(data=vardf, color='#dd1c77', alpha=0.8)
    sb.scatterplot(data=vardf, color='#df65b0', s=mean.values )
    plt.tight_layout()
    plt.savefig('../../output/manuscript/var_robustness_legend.png', dpi=300, transparent=True)
    plt.show()

    plt.figure(figsize=(9,2.5))
    sb.kdeplot(data=melt, x='value', hue='variable', fill=True, palette='RdPu')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/kde_robustness_legend.png', dpi=300, transparent=True)
    plt.show()

    plt.figure(figsize=(6,2.5))
    sb.kdeplot(data=melt, x='value', hue='variable', fill=True, palette='RdPu',legend=False)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/kde_robustness.png', dpi=300, transparent=True)
    plt.show()


    escapedf.to_csv('../../output/cocktail/cocktail_robustness.csv')




def graph_cocktail_distance():

    #calculate_cocktail_distance()
    #exit()

    figname = 'graph_cocktail_distance'

    ab_class_dict, ab_gene_dict,ab_label_dict = utils.get_antibody_properties()

    ''' Number of places that are significantly different between numabs ''' 
    allstats = pd.read_csv('../../output/cocktail/cocktail_differences.csv')
    allcocktailfit = pd.read_csv('../../output/cocktail/cocktail_virus_fitness.csv')
    neut = pd.read_csv('../../output/cocktail/cocktail_neut_fitness.csv')
    ace2 ='single_mut_effects_cleaned_with_predictors.csv'
    ace2df  = pd.read_csv('../../output/estimator/'+ace2)[['bind_avg','mutation']]


    ''' Graph ace2 histogram ''' 

    plt.figure(figsize=(2,2))
    sb.set(context='paper', font_scale=1.2, style='ticks')
    sb.histplot(data=ace2df, x="bind_avg", color='Red', alpha=0.7)
    plt.axvline(x=-0.5, lw=0.5, color='grey')
    plt.title('ACE2/RBD')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/ace2.png', dpi=300, transparent=True)

    abfit = utils.get_antibody_fitness_landscape().merge(ace2df, how='inner', left_on='mut', right_on='mutation')


    cocktails = ['C13','C15','C16']

    abcolors = fmt.CLASS_AB_COLOR_DICT

    for c in cocktails: 

        neuts = neut.loc[neut.cocktail ==c]
        cols = list(neuts.antibody) + ['mut', 'bind_avg']
        ss = abfit[cols].set_index('mut')

    
        abcolors= fmt.CLASS_AB_COLOR_DICT

        range_ab = [-3,4]
        range_ace2 = [-6,1.5]

        '''Repeat of KDE change''' 
        for col in neuts.antibody.values:
            
            abf = ss[[col, 'bind_avg']]
            mergedvariants = abf.loc[abf.index.isin(fmt.ALL_VARIANTS)]
            plt.figure(figsize=(2.5,2.5))
            sb.set(context='paper', font_scale=1.2, style='ticks')
            sb.kdeplot(data=abf, x='bind_avg', y=col, fill=True, color=abcolors[col], alpha=0.8)
            ax =sb.scatterplot(data = mergedvariants,x= 'bind_avg', y=col ,color='red', marker='o', alpha=0.6)
            #sb.scatterplot(data=abf, x='bind_avg', y=col,color='gray', marker='o', alpha=0.6)
            plt.axhline( y=0, lw=0.5, color='gray')
            plt.axvline( x=0, lw=0.5, color='gray')
            plt.xlim(range_ace2)
            plt.ylim(range_ab)
            plt.title(str(col))
            plt.tight_layout()
            plt.savefig('../../output/manuscript/' + figname + c+ '_' + col + 'kde_locs.png', dpi=300, transparent=True)
            #plt.show()
            print(mergedvariants)

        classes = [ab_class_dict[ab] for ab in neuts.antibody]
        labels = [ab_label_dict[ab] for ab in neuts.antibody]

        palette = [abcolors[ab] for ab in neuts.antibody]
        neuts['labels'] = labels

        figsize =(2,3)
        if c =='C13': 
            figsize = (2,2.75)
        plt.figure(figsize=figsize)
        sb.barplot(data=neuts, x ='labels', y = 'pc_neutralized',palette=palette)
        plt.xticks(rotation=90)
        plt.ylim([0,100])
        plt.title(c)
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_' + c + '_bar.png', dpi=300, transparent=True)


    allcocktailfit['location'] = allcocktailfit.mut.str[1:-1]
    print(allstats)

    ''' For 3 antibody cocktails ''' 
    for n in allstats.numab.unique():
        cocktails = allstats.loc[allstats.numab ==n].variable.unique() 

        cocktailfit = allcocktailfit.loc[allcocktailfit['cocktail'].isin(cocktails)]


        ''' Merge with ACE2 binding ''' 

        merged = cocktailfit.merge(ace2df, how='inner', left_on='mut', right_on='mutation')
        merged = merged.drop('mutation', axis=1)

        print(merged.head())

        range_ab = [-4,1.5]
        range_ace2 = [-6,1.5]

        variantsall = pd.DataFrame()

        for cocktail in merged.cocktail.unique():
            mergedss = merged.loc[merged.cocktail == cocktail]
            mergedvariants = mergedss.loc[mergedss.mut.isin(fmt.ALL_VARIANTS)]

            variantsall = pd.concat([mergedvariants, variantsall])
            x = 'bind_avg'
            y = 'max_neutralization'
            plt.figure(figsize=(2.5,2.5))
            sb.set(context='paper', font_scale=1.2, style='ticks')
            sb.kdeplot(data=mergedss, x=x, y=y, fill=True, color='gray', alpha=0.6)
            ax =sb.scatterplot(data = mergedvariants,x=x, y=y,color='red', marker='o', alpha=0.6)
            plt.axhline( y=0, lw=0.5, color='gray')
            plt.axvline( x=0, lw=0.5, color='gray')
            plt.xlim(range_ace2)
            plt.ylim(range_ab)
            plt.title(cocktail)
            plt.tight_layout()
            plt.savefig('../../output/manuscript/' + figname + '_' + cocktail + 'kde_locs.png', dpi=300, transparent=True)

            ''' Merged variants''' 
            plt.figure(figsize=(3,2.5))
            sb.set(context='paper', font_scale=1.2, style='ticks')
            sb.barplot(data=mergedvariants, x='mut', y=y, fill=True, color='red', alpha=0.8)
            plt.xticks(rotation=90)
            plt.title(cocktail)
            plt.tight_layout()
            plt.savefig('../../output/manuscript/' + figname + '_' + cocktail + '_variants.png', dpi=300, transparent=True)
        

        '''Escape variants together ''' 
        
        ss = variantsall.loc[variantsall.cocktail.isin(['C13','C15','C16'])]
        ss = ss.set_index('mut')

        plt.figure(figsize=(3,2.5))
        sb.kdeplot(data=ss, x=x, y=y, fill=True, palette='RdPu', alpha=0.4, hue='cocktail', legend=False)
        plt.axvline( x=0, lw=0.5, color='gray')
        plt.axhline( y=0, lw=0.5, color='gray')
        plt.xticks(rotation=90)
        plt.ylim([-1.5,1.5])
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_' + cocktail + '_variants-bar.png', dpi=300, transparent=True)


        ss = variantsall.loc[variantsall.cocktail.isin(['C17'])]
        ss = ss.set_index('mut')

        plt.figure(figsize=(3,2.5))
        sb.kdeplot(data=ss, x=x, y=y, fill=True, palette='Greens', alpha=0.4, hue='cocktail', legend=False)
        plt.axvline( x=0, lw=0.5, color='gray')
        plt.axhline( y=0, lw=0.5, color='gray')
        plt.xticks(rotation=90)
        plt.ylim([-1.5,1.5])
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_' + cocktail + '_variants-bar-few.png', dpi=300, transparent=True)



def graph_cocktail_distance_grouped():

    #calculate_cocktail_distance()
    #exit()

    figname = 'graph_cocktail_distance'
    ''' Number of places that are significantly different between numabs ''' 
    allstats = pd.read_csv('../../output/cocktail/cocktail_differences.csv')
    allcocktailfit = pd.read_csv('../../output/cocktail/cocktail_virus_fitness.csv')

    allcocktailfit['location'] = allcocktailfit.mut.str[1:-1]
    print(allstats)

    ''' For 3 antibody cocktails ''' 
    for n in allstats.numab.unique():
        cocktails = allstats.loc[allstats.numab ==n].variable.unique() 

        cocktailfit = allcocktailfit.loc[allcocktailfit['cocktail'].isin(cocktails)]
        cocktailfit = cocktailfit.loc[cocktailfit['max_neutralization']>fmt.MAX_AB_DDG_BIND]
            
        ''' Grouped ''' 
        count = cocktailfit.groupby(['location','cocktail']).count().reset_index()        
        mean = cocktailfit.groupby(['location','cocktail']).mean().reset_index()
        merged= count.merge(mean,on=['location','cocktail'], how='inner')
        merged = merged[['location', 'cocktail', 'mut', 'max_neutralization_y']]

        # 4 color categories 

        merged = merged.sort_values('location')
        merged = merged.set_index('location')

        f,ax = plt.subplots(figsize=(7,3.75))
        colors = dict(zip(cocktails, sb.color_palette('Set3', len(cocktails))))

        figsize=(7,3)
        plt.figure(figsize=figsize)
        sb.set(context='paper', font_scale=1, style='ticks')
        area = merged.mut*10
        sb.scatterplot(data = merged,x='location', y = 'max_neutralization_y',hue='cocktail', palette='Paired',s=area,alpha=0.8)
        #sb.lineplot(data = merged,x='location', y = 'max_neutralization_y',hue='cocktail', palette='Paired',alpha=0.8)
        plt.axhline(y=0,lw=0.5,color='grey')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig('../../output/manuscript/' + figname + '_mix_' + str(n) +'.png', dpi = 300, transparent=True)
      
    exit()

    grouped = data.groupby(['location','cocktail','numabs']).sum().reset_index()
    grouped['uncovered'] = 20- grouped.iscovered.values 
    
    print(grouped.head())

    f,ax = plt.subplots(figsize=(4.5,3.75))
    sb.set(context='paper', font_scale=1)
    #sb.barplot(data = embedding, x= 'countabs', y = 'pc_virus_coverage',hue=hue, palette= hue_colors, edgecolor='gray')
    #sb.violinplot(data = grouped, x= 'cocktail', y = 'iscovered', edgecolor='gray',scale='width')
    #sb.scatterplot(data = grouped, x= 'iscovered', y = 'max_neutralization', edgecolor='gray')
    val = 'max_neutralization' 
    val = 'uncovered' 

    x= 'numabs'
    #x= 'cocktail'
    sb.violinplot(data = data, x= x, y = val, color='white', scale='width')
    sb.stripplot(data = data, x= x, y = val, edgecolor='gray', s=3,hue='cocktail')

    #sb.stripplot(data = grouped, x= 'location', y = val, edgecolor='gray', jitter=True, s=3)
    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname + 'absvcov.png', dpi = 300, transparent=True)
    plt.show()

    abs = data[['cocktail','numabs']].drop_duplicates().set_index('cocktail')
    pccov = data[['cocktail', 'pc_virus_coverage']].drop_duplicates().set_index('cocktail')
    maxneut = data[['cocktail', 'max_neutralization']].drop_duplicates().groupby('cocktail').mean().reset_index().set_index('cocktail')





def graph_UMAP_cocktail():

    figname = 'graph_UMAP_cocktail'

    data = pd.read_csv('../../output/cocktail/cocktail_virus_fitness.csv')
    abs = data[['cocktail','numabs']].drop_duplicates().set_index('cocktail')
    pccov = data[['cocktail', 'pc_virus_coverage']].drop_duplicates().set_index('cocktail')
    maxneut = data[['cocktail', 'max_neutralization']].drop_duplicates().groupby('cocktail').mean().reset_index().set_index('cocktail')

    
    print(abs)
    print(pccov)
    print(maxneut.values)

    val = 'max_neutralization'
    #val = 'iscovered'
    pivot = pd.pivot(data=data, index='mut', columns='cocktail', values=val)
    
    print(len(pivot)) 
    print(len(pivot) -pivot.sum(axis=0))

    '''
    sb.heatmap(data=pivot,)
    plt.show()
    '''

    alg = 'UMAP'
    embedding = ma.learn_antibody_distance(pivot, alg)
    embedding['cocktail'] = abs.index
    embedding['countabs'] = abs.values
    embedding['pc_virus_coverage'] = pccov.values/pccov.values.max()
    embedding['maxneut'] = maxneut.values

    colors= 'Set3'
    hue = 'pc_virus_coverage'
    hue_colors = sb.color_palette('PuBu', len(pccov.values))     

    area = embedding.countabs**2
    f,ax = plt.subplots(figsize=(3.75,3.75))
    sb.set(context='paper', font_scale=1)
    #sb.barplot(data = embedding, x= 'countabs', y = 'pc_virus_coverage',hue=hue, palette= hue_colors, edgecolor='gray')
    sb.scatterplot(data = embedding, x= 'maxneut', y = 'pc_virus_coverage',hue=hue,s=area, palette= hue_colors, edgecolor='gray')

    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname + 'absvcov.png', dpi = 300, transparent=True)
    plt.show()


    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(3.75,3.75))
    sb.set(context='paper', font_scale=1)
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue, palette= hue_colors, legend=False, s=area, edgecolor='gray')

    for j, lab in enumerate(embedding.cocktail):
        ax.annotate(lab, (embedding.dim1[j]-0, embedding.dim2[j]+ 0.3))                

    plt.title(alg)
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname + val + 'annotation.png', dpi = 300, transparent=True)
    plt.show()

    sb.set(context='paper', style='white', font_scale=1.2)
    f,ax = plt.subplots(figsize=(3.75,3.75))
    sb.set(context='paper', font_scale=1)
    sb.scatterplot(data = embedding, x= 'dim1', y = 'dim2',hue=hue, palette= hue_colors, legend=False, s=area, edgecolor='gray')
    plt.title(alg)
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.tight_layout()
    plt.savefig('../../output/manuscript/' + figname + val + '.png', dpi = 300, transparent=True)
    plt.show()



def graph_cocktail_mixes_data(data):

    figname = 'optimal_cocktails'
    abs = data.columns[0:-8] 

    print(abs)
    ab_class_dict, ab_gene_dict,ab_label_dict = utils.get_antibody_properties()
    ablabels = [ab_label_dict[ab] for ab in abs]
    abclasses = [ab_class_dict[ab] for ab in abs]
    replace_dict = dict(zip(abs,ablabels))

    data = data.rename(columns = replace_dict)

    data = data.sort_values('num_abs', ascending=False)
    data = data.set_index('cocktail')

    covered = data.pop('inf_cov')
    numabs = data.pop('num_abs')
    coverage = data.pop('cov')
    isopt = data.pop('isopt')
    opt_cov = data.pop('opt_cov')
    neut_cov = data.pop('neut_cov')

    orreds = sb.color_palette('OrRd_r', len(covered)+5) 
    covered_color = dict(zip(covered.unique(), orreds))
    row_color1 = covered.map(covered_color)

    rdpu = sb.color_palette('Purples_r', 7) 
    numabs_color = dict(zip(numabs.unique(), rdpu))
    row_color2 = numabs.map(numabs_color)

    blues = sb.color_palette('Blues_r', len(coverage)) 
    coverage_color = dict(zip(coverage.unique(), blues))
    row_color3 = coverage.map(coverage_color)
    row_colors = [row_color3, row_color1, row_color2]

    sb.set(context='paper', style='white')
    g = sb.clustermap(data=data,cmap=fmt.INCOCKTAIL, linewidths= 0.5, linecolor='white', row_cluster=False, col_cluster=False, figsize=(4.5,5.5), row_colors=row_colors)
    plt.tight_layout() 
    plt.savefig('../../output/figs/' + figname + '.png', dpi=300, transparent=True)
    plt.show()


def graph_cocktail_mixes():

    figname = 'graph_cocktail_mixes'
    data = pd.read_csv('../../output/cocktail/cocktail_incidence_matrix.csv')

    print(data)
    abs = data.columns[0:16]
    ab_class_dict, ab_gene_dict,ab_label_dict = utils.get_antibody_properties()
    ablabels = [ab_label_dict[ab] for ab in abs]
    abclasses = [ab_class_dict[ab] for ab in abs]
    replace_dict = dict(zip(abs,ablabels))

    data = data.rename(columns = replace_dict)
    data = data.set_index('Cocktail')

    covered = data.pop('pc_virus_coverage')
    data = data.drop('mix',axis=1)

    purples = sb.color_palette('PuBu_r', len(covered)+5) 
    covered_color = dict(zip(covered.unique(), purples))
    row_colors = covered.map(covered_color)

    data = data.drop_duplicates()

    sb.set(context='paper', style='white')
    g = sb.clustermap(data=data,cmap=fmt.INCOCKTAIL, linewidths= 0.5, linecolor='white', row_cluster=False, col_cluster=False, figsize=(4,5), row_colors=row_colors)
    plt.tight_layout() 
    plt.savefig('../../output/manuscript/' + figname + '.png', dpi=300, transparent=True)
    plt.show()
    exit()

    cols = list(data.columns[0:16]) + list(['pc_virus_coverage', 'Number Abs','sum_neutralization'])
    data = data[cols]

    grouped = data.groupby('Cocktail').mean().reset_index()
    grouped['sum_neut_norm'] = grouped['sum_neutralization']/grouped['sum_neutralization'].min()

    data = data.reset_index()
    data = data.drop_duplicates()

    numabs = data.pop('Number Abs')
    sum_neut = data.pop('sum_neutralization')

    cols = grouped.columns[1:20]
    grouped = grouped[cols]
    grouped = grouped.sort_values('pc_virus_coverage', ascending=False)
    grouped['ylabel'] = grouped.pc_virus_coverage.astype(int).astype(str) + '%'
    grouped = grouped.set_index('ylabel')

    grouped = grouped.drop('pc_virus_coverage', axis=1)

    area = grouped.sum_neut_norm.values**4*300




def graph_optimization_results(cocktail_files, gamma1_range, gamma2_range):

    dfall = pd.DataFrame()
    ''' Get all files, calculate number abs used gamma range '''
    for f in cocktail_files: 
        df = pd.read_csv(f) 
        df = df.set_index('antibody') 
        dfsum = df.sum(axis=0)
        dfall = pd.concat([dfall, dfsum])

    rename_cols = dict(zip ([0],['numabs']))
    dfall = dfall.rename(columns=rename_cols)
    dfall = dfall.reset_index()

    cols = dfall['index'].str.split('_', expand=True)
    coverage = cols.iloc[:,0].astype(float)*100.
    dfall['coverage'] = np.round(coverage,0).astype(str) 
    dfall['gamma1'] = np.round(cols.iloc[:,1].astype(float),0)
    dfall['gamma2'] = np.round(cols.iloc[:,2].astype(float),0)

    dfall = dfall.set_index('index')

    ''' Show number of antibodies changes as a function of gamma1 '''
    figsize = (4,4)
    plt.figure(figsize= figsize)
    sb.set(context='paper', style='ticks', font_scale=1.3)
    ax1= sb.scatterplot(data=dfall, x= 'gamma1', y='numabs', hue='coverage', alpha = 0.4, marker='o', legend=False)
    ax2 = sb.lineplot(data=dfall, x= 'gamma1', y='numabs', hue='coverage', alpha = 0.9, legend=False)
    #plt.legend(ncol=3)
    plt.xticks(rotation = 90)
    plt.xlabel('$\gamma_1$')
    plt.ylabel('Number Ab in mix')
    plt.tight_layout() 
    plt.savefig('../../output/manuscript/combinations_gamma1.png', dpi=300)
    plt.show()    


def concatenate_cocktail_files(cocktail_files): 

    dfall = pd.DataFrame()
    for f in cocktail_files: 
        df = pd.read_csv(f)
        df = df.set_index('antibody')
        dfall = pd.concat([dfall, df], axis=1) 

    return dfall

''' Get every mix solved by coverage, For every mix in coverage, what is maximum coverage achieved '''

def get_antibody_mix_properties(cocktail_files):

    dfall = concatenate_cocktail_files(cocktail_files) 
    dft = dfall.transpose()
    dft = dft.drop_duplicates()

    tmp= ''
    for col in dft.columns[1:]: 
        tmp = tmp + dft[col].astype(int).astype(str)

    tmp = pd.DataFrame(tmp)
    tmp = tmp.drop_duplicates()
    indices = tmp.index

    mixes = dft.loc[dft.index.isin(indices)]

    cnames = ['C'+ str(i+1) for i in range(len(mixes.index))]
    mixes['Cocktail'] = cnames 
    mixes= mixes.reset_index()
    mixes = mixes.drop('index', axis=1)
    mixes= mixes.set_index('Cocktail')
    
    ''' Calculate maximum coverage for antibody mix ''' 
    data = utils.get_antibody_fitness_landscape()
    v = ['location', 'mut']

    neutall, mergedall = pd.DataFrame(), pd.DataFrame()

    for rown in range(len(mixes)):
        row = pd.DataFrame(mixes.iloc[rown,:])
        colname = row.columns[0]

        rowtf = row == True
        abs= rowtf.loc[rowtf[colname]==True].index.values
        cols = v+ list(abs)
        datamix = data[cols]

        neutv,countv= [],[]
        neut = pd.DataFrame()
        for ab in abs:
            totalmut = len(datamix[ab])
            totalmut = 2360
            d = datamix[ab].loc[datamix[ab]<=0]
            neutv.append(np.round(d.sum()))
            countv.append(np.round(d.count()/totalmut*100))
        
        neut['antibody'] = abs
        neut['sum_neutralization'] = neutv 
        neut['pc_neutralized'] = countv 
        neut['cocktail'] = colname

        #cols = pd.Series(neut.mix).str.split('_', expand=True).rename(columns={0:'cov', 1:'gamma1', 2:'gamma2'})
        #neut = pd.concat([neut,cols], axis=1)

        ''' Construct fitness matrix '''
        datamix['max_neutralization'] = datamix[datamix.columns[2:]].min(axis=1)

        ''' Merge with ACE2 binding ''' 
        ace2 ='single_mut_effects_cleaned_with_predictors.csv'
        ace2df  = pd.read_csv('../../output/estimator/'+ace2)[['bind_avg','mutation']]
        merged = datamix.merge(ace2df, how='inner', left_on='mut', right_on='mutation')
        merged = merged.drop('mutation', axis=1)

        ''' Incidence matrix ''' 
        incidence = merged[merged.columns[2:-2]] <= 0
        merged['iscovered'] = incidence.sum(axis=1)>0
        merged['cocktail'] = colname
        merged['numabs'] = len(abs)

        total_cov = merged.iscovered.sum()/len(merged.iscovered)*100
        neut['pc_virus_coverage'] = total_cov
        neutall = pd.concat([neutall, neut])       

        merged['pc_virus_coverage'] = total_cov
        cols = ['mut','max_neutralization', 'iscovered', 'cocktail', 'numabs','pc_virus_coverage']
        merged = merged[cols]
        mergedall = pd.concat([merged, mergedall])
        
    mergedall.to_csv('../../output/cocktail/cocktail_virus_fitness.csv', index=False)
    neutall.to_csv('../../output/cocktail/cocktail_neut_fitness.csv', index=False)

    mixes = mixes.reset_index()    
    merged = mixes.merge(neutall[['pc_virus_coverage','cocktail']], how='left', left_on='Cocktail', right_on='cocktail').drop('cocktail',axis=1).drop_duplicates()
    merged.to_csv('../../output/cocktail/cocktail_incidence_matrix.csv', index=False)


 

def graph_simulations_old(cocktail_file, gamma1, gamma2):

    ''' Graph antibody cocktail simulations '''
 
    allres = pd.read_csv(cocktail_file)
    allres = allres.replace({'C105_TH28I_YH58F':'C105$^{T28I/Y58F}$'})
    melted = allres.melt(id_vars='antibody', value_vars= allres.columns[1:])
    allres= allres.set_index('antibody')
    allres = allres.drop('C002-S2', axis=0)
    #allres = allres.drop('0.5', axis=1)

    splitcol = pd.DataFrame(melted['variable'].str.split(pat= '_', expand=True))
    splitcol = splitcol.rename(columns={0:'coverage', 1:'gamma1', 2:'gamma2'})
    melted = pd.concat([melted, splitcol], axis = 1)
    melted = melted.set_index('antibody')
    melted = melted.rename(columns={ 'variable':'params', 'value':'numabs'})                           
    grouped = melted.groupby(['params', 'coverage','gamma1']).sum().reset_index()
    grouped['gamma1'] = grouped.gamma1.astype(float)
    grouped = grouped.loc[grouped.gamma1 <= 1250]
    print(grouped.head())

    sorted = grouped.sort_values('gamma1')
    sorted['floatcov']  = sorted.coverage.astype(float)
    print(sorted.head())

    maxes = []
    abs = range(1,10)
    for ab in abs:
        maxes.append(sorted.loc[sorted.numabs == ab].floatcov.max())
        ss = sorted.loc[sorted.numabs == ab]
        print(ss)
        print('=====')
        print(allres.columns)
        allss  = allres[ss.params.unique()]
        print(allss.head())

        figsize = (5.5,5.5)
        cmap = ['#fde0dd','#c51b8a']
        plt.figure(figsize= figsize)
        sb.set(context='paper', style='white')
        sb.heatmap(data=allss ,cmap=cmap,square=True, linecolor='white', linewidths=0.5, cbar=None, alpha=0.7)
        plt.xticks(rotation = 90)
        plt.tight_layout() 
        plt.savefig('../../output/figs/combinations_' +str(ab) + '.png', dpi=300)
        plt.show()    
        
    maxdf = pd.DataFrame(maxes, columns={'max'})
    print(maxdf)


    ''' Show number of antibodies changes as a function of gamma1 '''
    figsize = (5.5,2.5)
    plt.figure(figsize= figsize)
    sb.set(context='paper', style='ticks')
    ax1= sb.barplot(data=maxdf, x=maxdf.index,y= maxdf['max'])
    plt.xticks(rotation = 90)
    plt.ylim([0,1])
    plt.ylabel('Max. input coverage')
    plt.xlabel('Count Ab')
    plt.tight_layout() 
    plt.savefig('../../output/figs/combinations_gamma1.png', dpi=300)
    plt.show()        


    ''' Show groups of antibodies, how they are picked, and max coverage, groups of three, four, three, two, one ... '''

    numabs = allres.sum(axis=0)
    allres = allres.transpose()
    allres['sum'] =allres.sum(axis=1)

    splitcol = allres.index.str.split(pat= '_', expand=True)

    splitcol = splitcol.rename(columns={0:'coverage', 1:'gamma1', 2:'gamma2'})
    allres = pd.concat([allres, splitcol], axis = 1)

    cols = allres.columns.astype(float)*100    
    cols = cols.astype(int).astype(str) +'%'
    #cols = allres.columns
    
    dictcols = dict(zip(allres.columns, cols))
    allres = allres.rename(columns=dictcols)

    allres = allres.transpose()

    figsize = (2.5,5.5)
    cmap = ['#fde0dd','#c51b8a']
    plt.figure(figsize= figsize)
    sb.set(context='paper', style='white')
    sb.heatmap(data=allres ,cmap=cmap, square=True, linecolor='white', linewidths=0.5, cbar=None, alpha=0.7)
    plt.xticks(rotation = 90)
    plt.tight_layout() 
    plt.savefig('../../output/figs/combinations.png', dpi=300)
    plt.show()    


def analyze_cocktail(antibody_fitness,virus_fitness, cocktail, cocktailname, graph=True):

    ''' Analyze cocktail in terms of its virus sensitvity 
        Graph virus location sensitivities 
    '''
    
    ''' Bar graph with the antibodies that are chosen virus coverage '''
    dfc = antibody_fitness[cocktail]
    dfc = dfc.reset_index()
    
    ''' Find min in rows and then pick those minimum positive ddgs => best chance at coverage '''
    dfc['min_ddg'] = dfc.iloc[:,1:].min(axis=1)
    dfc['ispos'] = dfc.min_ddg>0 
    
    cov = float(cocktailname)
    allmuts = len(dfc.mut.unique())
    uncovmuts = len(dfc.loc[dfc.ispos].mut.unique())
    covmuts = allmuts - uncovmuts 
    descov = allmuts*(cov)

    dfc['location'] = dfc.mut.str[1:-1]
    
    ''' All positive ispos ddgs are uncovered''' 
    dfu = dfc.loc[dfc.ispos]
    dfu['ddgn'] = dfu.min_ddg/dfu.min_ddg.max()    
    
    merged = vdf.merge(dfu, on='mut', how='inner')
    merged['minACE2']= -merged.ACE2
    merged['cocktail'] = cocktailname 

    ''' Scatterplot for each antibody '''
    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='white')
    ax = sb.scatterplot(data=merged, y='ddgn', x = 'minACE2', alpha=0.8,hue ='location', palette='Set3')
    plt.axhline(y=0, color='grey',xmin=-1, xmax=1, lw=0.5)
    plt.axvline(x=0, color='grey',ymin=-1, ymax=1, lw=0.5)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.ylabel('ddg(Ab/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_sensitivity' + cocktailname + '.png', dpi=300)

    ''' Do histogram/KDE plot '''
    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='white')
    ax = sb.histplot(data=merged,  x = 'minACE2', alpha=0.5)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_kde' + cocktailname + '.png', dpi=300)

    ''' Violin plot of all sensitive locations '''
    plt.figure(figsize=(3.5,1.5))
    sb.set(context='paper', style='ticks')
    sb.violinplot(data=merged, y='minACE2', x = 'location', color='white',scale='width', inner='point')
    ax = sb.stripplot(data=merged, y='minACE2', x = 'location', alpha=0.5, color='Red', jitter=True)
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_infectionsensitivity_' + cocktailname +'.png', dpi=300)

    if graph: 
        plt.show()
        
    return merged 


def compare_cocktails(antibody_fitness_file, virus_fitness_file, cocktail_file, cocktail_1, cocktail_2):

    ''' Compare cocktail sensitivities '''

    allres = pd.read_csv(cocktail_file)

    ''' Pick ones to compare ''' 

    allres = allres[['antibody',cocktail_1, cocktail_2]]

    c1 = allres.loc[allres[cocktail_1] == 1].antibody.values 
    c2 = allres.loc[allres[cocktail_2] == 1].antibody.values 

    abdf, vdf = import_fitness_landscapes(antibody_fitness_file, virus_fitness_file)

    c1data = analyze_cocktail(abdf, vdf,c1, cocktail_1, graph=False)
    c2data = analyze_cocktail(abdf, vdf,c2, cocktail_2, graph=False)

    merged = pd.concat([c1data, c2data], axis=0)

    colors = ['#de77ae', '#7fbc41']
    hue_colors = ['0.85','0.6']

    ''' Histogram/KDE plot '''
    plt.figure(figsize=(1.75,1.75))
    sb.set(context='paper', style='ticks')
    ax = sb.histplot(data=merged,  x = 'minACE2', alpha=1, hue='cocktail', palette= colors, hue_order=hue_colors)
    ax.legend().remove()
    plt.xlabel('ddg(ACE2/RBD)')
    plt.xticks(rotation=0)
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_kde.png', dpi=300)
    
    print(merged.head())
    sorted= merged.sort_values('location')

    print(sorted.head())

    plt.figure(figsize=(3.5,1.5))
    sb.set(context='paper', style='white')
    ax =sb.violinplot(data=sorted, y='minACE2', x = 'location', color='white',scale='width', inner='point', hue='cocktail', hue_order=hue_colors,palette=colors, legend=False, lw=0.5,alpha=0.7)
    #ax = sb.scatterplot(data=sorted, y='minACE2', x = 'location', alpha=0.5, palette=colors, hue='cocktail')
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_sensitivity.png', dpi=300)

    ''' Group the cocktails '''
    cmean = sorted.groupby(['location','cocktail']).mean().reset_index()[['cocktail','location','minACE2','ddgn']]
    cmean = cmean.loc[cmean.minACE2 <=0]
    sorted = sorted.loc[sorted.minACE2 <= 0]

    print(cmean.location.unique())
    plt.figure(figsize=(1.5,2.5))
    sb.set(context='paper', style='white')
    ax = sb.scatterplot(data=sorted, x='minACE2', y = 'location', alpha=0.8, palette=colors, hue='cocktail',hue_order=hue_colors, s=30, marker='s', size='min_ddg')
    #ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    #plt.xticks(rotation=90)
    #plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_sensitivity.png', dpi=300)
    plt.show()

    print(cmean.head())


    plt.figure(figsize=(2.5,2.5))
    sb.set(context='paper', style='ticks')
    #ax =sb.violinplot(data=sorted, y='minACE2', x = 'location', color='white',scale='width', inner='point', hue='cocktail', palette='Set2', legend=False, lw=0.5,alpha=0.7)
    ax = sb.kdeplot(data=sorted, x='minACE2', alpha=0.5, palette=colors, hue='cocktail')
    ax.legend().remove()
    plt.axhline(y=0, color='grey', lw=0.5)
    plt.xticks(rotation=90)
    plt.ylabel('ddg(ACE2)')
    plt.tight_layout() 
    plt.savefig('../../output/figs/combination_scatter.png', dpi=300)


def run_cocktail(args):
    parser = setup_parser()
    namespace = parser.parse_args()
    create_output_structure(namespace.output_dir)
    
    a_file = f'{namespace.filepath}/{namespace.filename}.csv'
    v_file = f'{namespace.vfilepath}/{namespace.vfilename}.csv'
    df, v_df = import_fitness_landscapes(a_file, v_file)


    #print(df.index == v_df.index)

    #print (df.head())
    #fitness_matrix = df.values
    #virus_fitness_matrix = v_df.values
    results = compute_cocktail(df, v_df, namespace.p_unmanaged, namespace.gamma1, namespace.gamma2)
    
    if results is None:
        warnings.warn('Maximum proportion of not managed virus is too low (no combination of antibodies can satisfy this constraint) so optimization not solved.')
    
    else:
        chosen_antibodies = df.columns[results.astype(bool)]
        with open(namespace.output_dir+'output/optimize/optimal_antibodies.txt', 'w') as f:
            f.write(','.join(chosen_antibodies))
 
    combinepath = f'{namespace.filename}_{namespace.vfilename}'
    with open(namespace.output_dir+'output/cocktail/'+combinepath+'_opt_results.txt', 'w') as writer:
        for choice in choices:
            writer.write(choice+'\n')




