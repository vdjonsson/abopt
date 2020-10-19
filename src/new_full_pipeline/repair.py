import energy as e
import foldx as foldx
import pandas as pd
import structure
import os

def repair (pdb_dirs, pdb_list, out_dirs):

    i = 0
    for pdb in pdb_list:
        pdb_dir = pdb_dirs[i]
        out_dir = out_dirs[i]
        if os.path.isdir(out_dir) == False:
            os.mkdir(out_dir)

        foldx.run_repair_model(pdb, pdb_dir, out_dir)
        i= i+1 
