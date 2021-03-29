import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
import plotutils as pu 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.stats as stats 
from sklearn import preprocessing


filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'


def config():

    ''' Read ab configuration file '''
    abdf = pd.read_table('../../manuscript/antibody_list.txt', sep=',') 
    
    abdf= abdf.loc[~abdf.antibody.isin(['CR3022','C104', 'C121-S2','S309'])]
    ab_names = abdf.antibody.astype(str).values
    pdb_names = abdf.pdb.astype(str).values
    rbd_chains = abdf.rbdchain.astype(str).values

    pdb_rbd =dict(zip(pdb_names, rbd_chains))
    ab_pdb = dict(zip(ab_names, pdb_names))

    pdb_files = dict(zip(pdb_names, [pdb +'.pdb' for pdb in pdb_names]))
    pdb_dirs = dict(zip(pdb_names,[ '../../data/pdb/' + ab +'/unrepaired/' for ab in ab_names]))
    repair_dirs = dict(zip(pdb_names, [ '../../output/repair/' + ab +'/' for ab in ab_names]))
    remove_dirs = dict(zip(pdb_names,  [ '../../output/remove/' + ab +'/' for ab in ab_names]))
    constrain_dirs =dict(zip(pdb_names, [ '../../output/constrain/' + ab +'/' for ab in ab_names]))
    scan_dirs = dict(zip(pdb_names,[ '../../output/scan/' + ab +'/' for ab in ab_names]))
    energy_dirs = dict(zip(pdb_names,[ '../../output/energy/' + ab +'/' for ab in ab_names]))
    design_dirs = dict(zip(pdb_names,[ '../../output/design/' + ab +'/' for ab in ab_names]))
    mutate_dirs = dict(zip(pdb_names,[ '../../output/mutate/' + ab +'/' for ab in ab_names]))

    return ab_names, pdb_names, pdb_rbd, ab_pdb, scan_dirs, energy_dirs

def ddgbind():

    ab_names, pdb_names, pdb_rbd, ab_pdb, scan_dirs, energy_dirs = config()

    ''' Calculate and merge Gibbs energy of binding '''     
    for ab_name in ab_names: 

        pdb_name = ab_pdb[ab_name]

        ' Set up files '
        pdb = ab_pdb[ab_name] + '_Repair'
        pdb_less = ab_pdb[ab_name] + '_Repair_less_ab_Repair'

        ddg = energy.calculate_ddg_bind(ab_name, pdb, pdb_less, scantype='virus', indir= scan_dirs[pdb_name], outdir= energy_dirs[pdb_name])


def merge():        

    ab_names, pdb_names, pdb_rbd, ab_pdb, scan_dirs, energy_dirs = config()
    merge_dir  = '../../output/merge/'
    merged = pd.DataFrame()

    i = 0 
    ''' Merge data for all antibodies ''' 
    for ab_name in ab_names: 

        pdb_name = ab_pdb[ab_name]
        energy_dir = energy_dirs[pdb_name]
        df = pd.read_csv(energy_dir + 'ddg_bind_' + ab_name + '_virus_scanning.csv')

        df = df[['mut', 'ddg_bind']]
        df = df.rename(columns = {'ddg_bind':ab_name})

        if i > 0:
            merged = df.merge(merged, how= 'outer')
        else:
            merged = df

        print (merged)
        i = i+1

    merged = merged.loc[~merged.mut.str[-1].isin(['e', 'o'])]
    merged['location'] = merged.mut.str[1:-1].astype(int)
    merged= merged.loc[merged['location'] >=400]
    merged.to_csv(merge_dir + 'rbd_ab_fitness.csv', index=None)

    ''' Output the Ace2/RBD landscape '''     
    rbd = pd.read_csv('../../data/estimator/single_mut_effects_cleaned.csv')
    
    rbd = rbd[['mutation','bind_avg']]
    rbdneg  = rbd.loc[rbd.bind_avg <=0]
    rbdpos  = rbd.loc[rbd.bind_avg >0]
    maxval = rbdpos.bind_avg.max()
    minval = rbdneg.bind_avg.min()

    rbdneg['bind_norm'] = -rbdneg.bind_avg/minval
    rbdpos['bind_norm'] = rbdpos.bind_avg/maxval
 
    rbd_merged = pd.concat([rbdneg, rbdpos],axis=0)
    rbd_merged = rbd_merged.rename(columns={'bind_norm':'ACE2'})

    print(rbd_merged.head())
    rbd_merged.to_csv(merge_dir + 'rbd_ace2_fitness.csv', index=None)

    ''' Merge fitness landscapes '''
    merged_fit = merged.merge(rbd_merged, how='inner', left_on='mut', right_on='mutation')
    merged_fit = merged_fit.drop(['mutation', 'bind_avg', 'location'], axis=1)
    ab_fit =  merged_fit.iloc[:,0:-1]
    rbd_fit =  merged_fit[['mut', 'ACE2']]

    ab_fit.to_csv(merge_dir + 'rbd_ab_fitness.csv', index=None)
    rbd_fit.to_csv(merge_dir + 'rbd_ace2_fitness.csv', index=None)
    merged_fit.to_csv(merge_dir + 'rbd_ace2_ab_fitness.csv', index=None)


ddgbind()    
merge()



