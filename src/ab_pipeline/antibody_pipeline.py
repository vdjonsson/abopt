import pandas as pd
import seaborn as sb 
import os 

data_path = '../../data/'

def read_estimator(filename, ab):
    df = pd.read_csv(data_path + filename)
    dfss = df.loc[df.antibody_id == ab ]
    return dfss 

def construct_positionscan_string (pdb_name, ab_pos, pdb_loc):

    merged = ab_pos.join(pdb_loc, on='fasta_location', how = 'left', lsuffix='_abpos')
    posscan = merged.pdb_wt_foldx.unique()
    pd.Series(posscan).to_csv(data_path + pdb_name + '_position_scan.csv', index=False)
    posscan_str =  ",".join(posscan)
    return posscan_str


# get positions where wildtype== True and coeff>=0, or wildtype==False, and coeff<=0 
def get_mutation_positions(ab):

    # get positions where wildtype== True and coeff>=0 
    mut1 = ab.loc[(ab.wild_type == True) & (ab.coefficient >=1e-10)]
    mut2 = ab.loc[(ab.wild_type == False) & (ab.coefficient <= -1e-10)] 
    abs = pd.concat([mut1, mut2])
    return abs 
    

def read_pdb_locations(ab_name):
    df = pd.read_csv(data_path + ab_name + '_locations.csv', dtype='object')
    df['pdb_wt_foldx'] = df.wt + df.pdb_chain + df.pdb_location + 'a'
    pdb = df.pdb.unique()[0]
    return pdb, df


def run_position_scan (ab_name, pdb_name, pos_scan):

    print(pdb_name)
    command = "foldx --command=PositionScan --pdb="+pdb_name+"_Repair.pdb "
    command = command + "--positions="+pos_scan +" --out-pdb=false"
    os.system(command)



f = 'nussenzweig_antibody_data_cleaned_with_alignments_mapped_back_sars_cov_2_ic50_ngml_coefficients.csv'
ab_name = 'C105'
ab = read_estimator(f,ab_name)

pdb_name, pdb_loc = read_pdb_locations(ab_name)
ab_pos = get_mutation_positions(ab)
pos_scan = construct_positionscan_string (pdb_name, ab_pos, pdb_loc)

print (pos_scan)

#run_position_scan (ab_name, pdb_name, pos_scan)
