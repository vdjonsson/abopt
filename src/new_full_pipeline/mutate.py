import energy as e
import foldx as foldx
import pandas as pd
import structure
import os

def mutate (pdb, mutations, pdb_dir, out_dir):

    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)

    foldx.create_individual_list(mutations,pdb,out_dir)
    foldx.run_build_model(pdb, 'individual_list.txt', pdb_dir, out_dir)
    foldx.rename_buildmodel_files(pdb[:-4], out_dir, './individual_list.txt')
