import pandas as pd 



def write(data, ismatrix=False):
    if ismatrix: 
        print(data.head())
    data.to_csv('../tmp/tmp.csv')
