import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import structure 
import energy     
import matplotlib.pyplot as plt 

filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'
file_ab_locations = '../../data/location/C105_locations.csv'


dpi = 300 

ab_names = ['B38', 'C105','CC121', 'CB6', 'COVA2-39','CV30']
pdb_names = ['7bz5', '6xcm', '6xc3', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'C', 'A', 'A', 'E']

ab_names = ['C105','B38', 'CB6', 'COVA2-39','CV30']
class_names= ['class I', 'class I', 'class I?', 'class II', 'class I?']
pdb_names = [ '6xcm', '7bz5', '7c01', '7jmp', '6xe1']
rbd_chains = ['A', 'C', 'A', 'A', 'E']


ab_names = ab_names[2:]
pdb_names = pdb_names[2:]
class_names = class_names[2:]
rbd_chains = rbd_chains[2:]

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




def combine_data():

    merged = pd.DataFrame()
    for ab in ab_names:
        title = ab + ' viral scanning ' + ab_class[ab]
        pdb = ab_pdb[ab]
        energy_dir = energy_dirs[pdb]

        filen ='ddgbind_' + ab + '_virus_scanning.txt'
        df = pd.read_table(energy_dir + filen, ',')
        df['antibody'] = ab
        df['mut_nochain'] = df.mut.str[0] + df.mut.str[2:]

        #pdb_file_name = pdb + '.pdb'
        #labeled_chains = structure.label_chains(pdb)
        #ep = structure.find_epitopes(pdb_dirs[pdb], pdb_file_name, labeled_chains, 4)    

        epitope = pd.read_csv('../../data/location/' + pdb + '_epitopes.csv')        

        proximity = epitope.number_virus.unique()

        # remove chain and mutation that end with e or o

        locations = pd.DataFrame(df.pdb_location.unique())

        print(locations)
        vlines = locations.loc[locations[0].isin(proximity)].index

        #sb.violinplot(data=df, x='pdb_location', y='ddg_bind')
        plt.figure(figsize=(12,2.5))
        sb.set(context='paper', style='ticks')
        sb.stripplot(data=df, x='pdb_location', y='ddg_bind', color='red', alpha=0.6)
        plt.vlines(vlines, ymin = df.ddg_bind.min(), ymax= df.ddg_bind.max(), color='#ffeda0', alpha=0.6, lw=2)
        plt.xticks(rotation=90)
        plt.xlabel('RBD location')
        plt.title(title)
        sb.despine()
        plt.tight_layout()
        plt.savefig(fig_dir + title + '.png', dpi =dpi)
        plt.show()
        
        merged = pd.concat([df, merged])


    merged.to_csv(merge_dir + 'ddgbind_virus_scanning.csv')
    print(merged)


#viral_scanning_binding()
combine_data()
