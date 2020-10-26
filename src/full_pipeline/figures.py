import pandas as pd 
import seaborn as sb 
import analysis as aa 
import antibody_pipeline as ap
#import utils as u 
import plotutils as pltab
import numpy as np           
import matplotlib.pyplot as plt 

out_fig ='../../output/figs/'
out_tmp = '../../output/tmp/' 
ab_names = ['CB6', 'C105', 'B38','CC12.1','ACE2']
pdb_names =['7c01','6XCM', '7bz5','6xc3', '6m0j']
has_mutations =[True, True, False, False]

ab_muts = dict(zip(ab_names, has_mutations))
ab_pdb = dict(zip(ab_names, pdb_names))


def parse_mutation(data): 

    data['mutation'] = data.mutation.values
    data['location'] = data.mutation.str[1:-1]
    data['mut'] = data.mutation.str[-1]

    return data 

def set_location_index(data):

    locations = data.location.unique()
    loc_dict = dict(zip(locations, range(len(locations))))
    location_index = [loc_dict[location] for location in data.location.values]
    data['location_index'] = location_index
    return data 

'''
def figure_1():  # plot fitness as violinplot with matrix data 

    fit = pd.read_csv(out_tmp + 'matrix.csv')    
    fit = parse_mutation(fit)

    # Plot ab correlation 
    abs = fit.columns[1:-2]

    fits = fit[abs]
    #u.write(fits, True)
    pltab.plot_correlation(data=fits, fignum='1a')
   
    ddgsm = pd.melt(ddgsall, id_vars=['location', 'mut', 'location_index'], value_vars=ddgsall.columns[0:4], var_name= 'ab',value_name='ddg')

    # get epitopes 
    epitopelines = get_epitope_lines(epitopename, epitopefile)
    epitopes = pd.read_csv('../../data/location/epitopes.csv')   
    epitopes['location'] = epitopes.epitope_location.astype(str)
    ddgsm = ddgsm.merge(epitopes, how = 'left', left_on='location', right_on='location')
    locs = pd.DataFrame(ddgsm.location.unique())
    rbdlines = ddgsm.loc[ddgsm.epitope_binding_molecule == 'ACE2'].location_index.unique()

    for ab in ab_names: 
        pltab.plot_ddg_stripplot(ddgsm,x='location', y='ddg', hue='ab', filtername='ab', filterval=ab, title='virus scan', epitopes=[rbdlines], epitopetype=['RBD'])    
        
    pltab.plot_ddg_stripplot(ddgsm,x='location', y='ddg', hue='ab',title='virus scan', epitopes=[rbdlines], epitopetype=['RBD'])    
'''
def figure_2(): # graph individual wt ab ddgs violin plots

    ab_names = ['C105','ACE2']
    ddgab = aa.combine_antibody_binding_energy_data(ab_names, scantype='virus')
    ddgab['location']= ddgab.location.astype(int)
    ddgab = ddgab.loc[(ddgab.location >= 400)]
    ddgab = ddgab.loc[(ddgab.location <= 520)]
    #u.write(ddgab)

    epitopes = pd.read_csv('../../data/location/epitopes.csv')
    epitopes['location'] = epitopes.epitope_location.astype(int)
    ddgab = ddgab.merge(epitopes, how = 'left', left_on='location', right_on='location')
    locs = pd.DataFrame(ddgab.location.unique())
    ddgab = set_location_index(ddgab)
    rbdlines = ddgab.loc[ddgab.epitope_binding_molecule == 'ACE2'].location_index.unique()

    for ab in ddgab.abname.unique():
        pltab.plot_ddg_stripplot(ddgab,x='location', y='ddg_bind', hue='abname', 
                                 filtername='abname', filterval=ab, title='virus scan',
                                 epitopes=[rbdlines], epitopetype=['RBD'])

    print(ddgab.columns)
    ddgab = ddgab[['mutation', 'ddg_bind','abname']]
    ddgab['ddg_bind_p'] = ddgab.ddg_bind +0.5
    pivot = ddgab.pivot(index = 'abname', columns='mutation', values='ddg_bind_p').transpose()



    sb.scatterplot(data=pivot, x='ACE2', y='C105')

    plt.show()
    print(pivot.head())
 
figure_2()
