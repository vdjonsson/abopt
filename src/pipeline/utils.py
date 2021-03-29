import pandas as pd 
from sklearn import preprocessing

aa_three = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE', 'LEU', 'LYS', 'MET','PHE', 'PRO','SER', 'THR','TRP', 'TYR','VAL','NA']
aa_single = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F', 'P', 'S', 'T','W', 'Y','V' ,'X']
aa_1to3_dic = dict(zip(aa_single, aa_three))
aa_3to1_dic = dict(zip(aa_three, aa_single))


def get_antibody_data(file):
    abdata = pd.read_table(file, sep=',')
    return abdata


def clean_filter_landscape(data): 

    locations = data.mut.str[1:-1].astype(int)
    data['location'] = data.mut.str[1:-1].astype(int)
    data = data.loc[data['location'] <520]
    ab_order = ['mut','location','C105','C105_TH28I_YH58F','B38','CB6','CV30','CC121','C002-S1', 'C002-S2','C119','C121-S1','C144','COVA2-39','C135','C110' ,'REGN10987','REGN10933']
    data = data[ab_order]
    return data 

def normalize_landscape(data, normalization = None):
    
    ''' Normalize antibodies '''
    ''' Scale between -1 and 1 for every antibody'''

    data = data.fillna(0)
    values = data.iloc[:,2:].values


    if normalization == None: 
        data_norm = data.iloc[:,1:]
        data_norm = data_norm.set_index('location')
        normalization = 'None'
    elif normalization in ['l1', 'l2','max']:
        normalizer = preprocessing.Normalizer(normalization)
        data_norm = pd.DataFrame(normalizer.transform(values), columns=data.columns[2:], index=data['location'])
    elif normalization == 'scale':
        ''' To do implement this '''
        tmp = pd.DataFrame()
        for ab in abs:
            s  = grouped.loc[grouped.antibody == ab]
            maxval  = s.binding.values.max()
            minval  = s.binding.values.min()
            sn = s.loc[s.binding<=0]
            sp = s.loc[s.binding>0]

            sn['binding'] = -sn.binding/minval 
            sp['binding'] = -sp.binding/maxval 
        
            tmp = pd.concat([sn, sp, tmp])    


    data_norm['mut'] = data.mut.values    

    data_norm.to_csv('../../output/merge/rbd_ab_fitness_normalized_' + normalization +'.csv')
    return data_norm
    

def get_antibody_fitness_landscape(normalization=None): 

    merge_dir = '../../output/merge/'
    data = pd.read_csv(merge_dir + 'rbd_ab_fitness_opt.csv')

    data_clean = clean_filter_landscape(data) 
    #data_opt = calculate_C105_optimized_landscape(data_clean)
    data_norm = normalize_landscape(data_clean, normalization)
    data_norm = data_norm.reset_index()

    return data_norm


def get_antibody_properties():

    abdf = pd.read_table('../../data/meta/antibody_list.txt', sep=',')
    ab_class = abdf['abclass'].values 
    ab_names = abdf['antibody'].values 
    ab_genes = abdf['vhgene'].values 
    ab_labels = abdf['label'].values 

    ab_class_dict = dict(zip(ab_names, ab_class))
    ab_gene_dict = dict(zip(ab_names, ab_genes))
    ab_label_dict = dict(zip(ab_names, ab_labels))

    return ab_class_dict, ab_gene_dict,ab_label_dict 
