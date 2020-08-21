import pandas as pd 
from biopandas.pdb import PandasPdb



data_path = '../../data/'
estimator_path = data_path + 'estimator/'
location_path = data_path + 'location/'
pdb_path = data_path + 'pdb/'
out_path = data_path + 'ddg/'

ab_name = 'CB6'
pdb_name = '7c01_Repair_AH60Y.pdb'

pdb = pb.read_pdb(pdb_path + ab_name +'/' +pdb_name)  

print('PDB Code: %s' % pdb.code)

