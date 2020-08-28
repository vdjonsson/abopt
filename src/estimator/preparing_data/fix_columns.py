import numpy as np
import pandas as pd

filepath = '../../data/'
filename = 'NeutSeqData_VH3-53_66_aligned'#'NeutSeqData_C002-215_cleaned_aligned'

correct_seq_col = 'sequence'
correct_hc_aligned = 'heavy_chain_aligned'
correct_lc_aligned = 'light_chain_aligned'
correct_hc = 'heavy_chain'
correct_lc = 'light_chain'
correct_id_col = 'antibody_id'

correct_cols = [correct_seq_col, correct_hc_aligned, correct_lc_aligned, correct_hc, correct_lc, correct_id_col]

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)

current_seq_col = 'sequence'#'sequence'
current_hc_aligned = 'VH or VHH_aligned'#'igh_vdj_aa_aligned'
current_lc_aligned = 'VL_aligned'#'igl_vj_aa_aligned'
current_hc = 'VH or VHH'#'igh_vdj_aa'
current_lc = 'VL'#'igl_vj_aa'
current_id_col = 'Name'#'antibody_id'

current_cols = [current_seq_col, current_hc_aligned, current_lc_aligned, current_hc, current_lc, current_id_col]

df.rename(columns = dict(zip(current_cols, correct_cols)), inplace=True)

df.to_csv(filepath+filename+'.csv', sep=',', header=True, index=False)
