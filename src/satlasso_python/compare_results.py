import pandas as pd
import numpy as np
import csv

def adjust_positions(a, offset):
    new_elements=[]
    for element in a:
        pos = int(element[1:len(element)])+offset
        new_elements.append(element[0]+str(pos))
    return new_elements

# def adjust_positions(dictionary, offset):
#     new_keys = []
#     for key,value in dictionary.items():
#         pos = int(key[1:len(key)])+offset
#         new_keys.append(key[0]+str(pos))
#     return new_keys

filepath = '../data/'
filename = 'single_mut_effects_cleaned'#'kyratsous_neutralization_data'
col = 'bind_avg'
df = pd.read_csv(filepath+filename+'_coefficients_'+col+'.csv', sep=',', header=0)
non_zero_coefficients = df.columns[df.iloc[0,:]!=0].values
non_zero_coefficients_values = df[non_zero_coefficients]
non_zero_coefficients_values.columns = adjust_positions(non_zero_coefficients_values.columns.values, 330-12)
non_zero_coefficients_values.to_csv(filepath+filename+'_aa_positions.csv', sep=',', header=True, index=False)

# df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0)

# dictionary = {}
# for col in df.columns[1:len(df.columns)]:
#     col_df = pd.read_csv(filepath+filename+'_coefficients_'+col+'.csv', sep=',', header=0)
#     non_zero_coefficients = col_df.columns[col_df.iloc[0,:]!=0].values
#     dictionary[col] = non_zero_coefficients

# aa_pos_array = []
# for key,val in dictionary.items():
#     aa_pos_array = aa_pos_array+val.tolist()

# aa_pos_set = set(aa_pos_array)
# count_dict = {}
# for item in aa_pos_set:
#     count_dict[item] = 0
#     for key,val in dictionary.items():
#         if item in val:
#             count_dict[item] = count_dict[item]+1

# new_keys = adjust_positions(count_dict, -19)
# count_dict = dict(zip(new_keys,list(count_dict.values())))
# file = open(filepath+filename+'_aa_positions.csv', 'w')
# writer = csv.writer(file)
# for key,value in count_dict.items():
#     writer.writerow([key,value])
# file.close()

# N354D, I472V most influential
