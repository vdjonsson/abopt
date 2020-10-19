import energy as e
import foldx as foldx
import pandas as pd
import structure
import os 

def scan (scantype, scanvalues, scanmolecule, antibody, pdblist, pdbdir, outdir):

    """ Mutational scanning
    """
    if os.path.isdir(outdir) == False:
        os.mkdir(outdir)

    if scantype == 'mutation':
        pdbname, pdb_loc = e.read_pdb_locations(file_location=file_ab_locations)
    elif scantype =='location':
        posscan_str =  ",".join(scanvalues)
    elif scantype == 'chain': # position scan entire chain
        posscan_str =  ","

    for pdb in pdblist:
        foldx.run_position_scan (pdb, scanmolecule,  posscan_str, pdbdir, outdir)
