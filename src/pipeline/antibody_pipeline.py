import energy as e 
import foldx as foldx
import pandas as pd
import structure 
import os 

def constrain(constraintype, constrainfile, antibody, cutoff , top, out_dir):

    """ Constrain data in constrainfile using cutoff 
    """
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
 
    if constraintype == 'estimator':
        estimator = e.read_estimator(constrainfile, antibody)
        constrained = e.constrain_estimator_features(estimator, cutoffmin =cutoff[0], cutoffmax=cutoff[1], topmuts = top, filter_name='chain', filter='H')
        constrained.to_csv(out_dir + antibody + '_estimator.csv', index=None)
    elif constraintype == 'energy':
        energies = pd.read_table(constrainfile, sep=',')
        constrained = e.constrain_energies(energies, cutoffmin= cutoff[0], cutoffmax= cutoff[1], topmuts=top)
        constrained.to_csv(out_dir + antibody + '_energies.csv', index=None)


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


def energy (antibody, pdb, pdb_less, scantype, energy_type, indir, outdir):

    """ Calculate energies 
    """ 

    if os.path.isdir(outdir) == False:
        os.mkdir(outdir)

    print(pdb)
    print(pdb_less)
    if energy_type =='ddgbind':
        ddg = e.calculate_ddg_bind(antibody,pdb, pdb_less, scantype=scantype, indir=indir, outdir=outdir)

    ddg.to_csv(outdir + energy_type + '_' + antibody + '_' + scantype + '_scanning.txt', index=None)



def design (designtype, designval, file_estimator, file_energies, out_dir):

    """ Returns places to mutate on antibody 
    """

    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)

    if designtype == 'antibody design':
        estimator = pd.read_csv(file_estimator)
        energies = pd.read_csv(file_energies)
        merged = energies.merge(estimator,how='inner', left_on='wt', right_on='wt_pdb')
        merged.to_csv(out_dir + designval + '_design.csv', index=None)
 

def mutate (pdb, mutations, pdb_dir, out_dir):

    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)

    foldx.create_individual_list(mutations,pdb,out_dir)
    foldx.run_build_model(pdb, 'individual_list.txt', pdb_dir, out_dir)
    foldx.rename_buildmodel_files(pdb[:-4], out_dir, './individual_list.txt')


def repair (pdb_dirs, pdb_list, out_dirs):

    i = 0 
    for pdb in pdb_list:
        pdb_dir = pdb_dirs[i]
        out_dir = out_dirs[i]
        if os.path.isdir(out_dir) == False:
            os.mkdir(out_dir)

        foldx.run_repair_model(pdb, pdb_dir, out_dir)
        i= i+1 

def remove(pdb_dirs, pdb_list, chains, chain_type, out_dirs):

    i = 0 
    for pdb in pdb_list:
        out_dir = out_dirs[i]
        pdb_dir = pdb_dirs[i]
        
        if os.path.isdir(out_dir) == False:
            os.mkdir(out_dir)
        pdb_less = structure.remove_chains(pdb_dir,pdb, chains, chain_type, out_dir)
        i = i+1
