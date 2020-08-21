import pandas as pd 
import seaborn as sb 
import antibody_analysis as aa 
import antibody_pipeline_vdj as ap
import utils as u 
import plot_utils as pltab
import numpy as np           
import matplotlib.pyplot as plt 

out_fig ='../../output/figs/'
out_tmp = '../../output/tmp/' 
ab_names = ['CB6', 'C105', 'B38','CC12.1']
pdb_names =['7c01','6XCM', '7bz5','6xc3']
has_mutations =[True, True, False, False]

ab_muts = dict(zip(ab_names, has_mutations))
ab_pdb = dict(zip(ab_names, pdb_names))


def parse_mutation(data): 

    data['mutation'] = data.mutation.values
    data['location'] = data.mutation.str[1:-1]
    data['mut'] = data.mutation.str[-1]

    # create a location index unique for the location
    locations = data.location.unique()
    loc_dict = dict(zip(locations, range(len(locations))))
    location_index = [loc_dict[location] for location in data.location.values]
    data['location_index'] = location_index
    return data 



def figure_1():  # plot fitness as violinplot with matrix data 

    fit = pd.read_csv(out_tmp + 'matrix.csv')    
    fit = parse_mutation(fit)

    # Plot ab correlation 
    abs = fit.columns[1:-2]

    fits = fit[abs]
    u.write(fits, True)
    pltab.plot_correlation(data=fits, fignum='1a')


    exit()
   
    ddgsm = pd.melt(ddgsall, id_vars=['location', 'mut', 'location_index'], value_vars=ddgsall.columns[0:4], var_name= 'ab',value_name='ddg')

    epitopes = pd.read_csv('../../data/location/epitopes.csv')   
    epitopes['location'] = epitopes.epitope_location.astype(str)
    ddgsm = ddgsm.merge(epitopes, how = 'left', left_on='location', right_on='location')
    locs = pd.DataFrame(ddgsm.location.unique())
    rbdlines = ddgsm.loc[ddgsm.epitope_binding_molecule == 'ACE2'].location_index.unique()

    for ab in ab_names: 
        pltab.plot_ddg_stripplot(ddgsm,x='location', y='ddg', hue='ab', filtername='ab', filterval=ab, title='virus scan', epitopes=[rbdlines], epitopetype=['RBD'])    
        
    pltab.plot_ddg_stripplot(ddgsm,x='location', y='ddg', hue='ab',title='virus scan', epitopes=[rbdlines], epitopetype=['RBD'])    


def figure_2():

    ddgs = aa.combine_antibody_binding_energy_data(ab_names, scantype='virus')
    u.write(ddgs)

figure_2()



