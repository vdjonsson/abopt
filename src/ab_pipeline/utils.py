import pandas as pd 



def write(data, filename = 'tmp', ismatrix=False):
    if ismatrix: 
        print(data.head())
    data.to_csv('../tmp/' + filename + '.csv')
