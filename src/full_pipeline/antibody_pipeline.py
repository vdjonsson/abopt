import energy as e 
import foldx as foldx
import pandas as pd
import structure 



def constrain(constraintype, constrainfile, antibody, cutoff , top):

    """ Constrain data in constrainfile using cutoff 
    """

    output_dir = '../../output/constrain/'
    if constraintype == 'estimator':
        estimator = e.read_estimator(constrainfile, antibody)
        constrained = e.constrain_estimator_features(estimator, cutoffmin =cutoff[0], cutoffmax=cutoff[1], topmuts = top, filter_name='chain', filter='H')
        constrained.to_csv(output_dir + antibody + '_estimator.csv', index=None)
    elif constraintype == 'energy':
        energies = pd.read_table(constrainfile, sep=',')
        constrained = e.constrain_energies(energies, cutoffmin= cutoff[0], cutoffmax= cutoff[1])
        constrained.to_csv(output_dir + antibody + '_energies.csv', index=None)


def scan (scantype, scanvalues, scanmolecule, antibody, pdbfile, pdbdir, outdir):

    """ Mutational scanning 
    """

    if scantype == 'mutation': 
        pdbname, pdb_loc = e.read_pdb_locations(file_location=file_ab_locations)
    elif scantype =='location':
        posscan_str =  ",".join(scanvalues)
    elif scantype == 'chain': # position scan entire chain
        posscan_str =  ","

    foldx.run_position_scan (pdbfile, scanmolecule,  posscan_str, pdbdir, outdir)


def energy (antibody, pdb, pdb_less, scantype, energy_type, indir, outdir):

    """ Calculate energies 
    """ 
    if energy_type =='ddgbind':
        ddg = e.calculate_ddg_bind(antibody,pdb, pdb_less, scantype='ab', indir=indir, outdir=outdir)

    ddg.to_csv(outdir + energy_type + '_' + antibody + '_' + scantype + '_scanning.txt', index=None)



def design (designtype, designval, file_estimator, file_energies):

    """ Returns places to mutate on antibody 
    """

    if designtype == 'antibody design':
        estimator = pd.read_csv(file_estimator)
        energies = pd.read_csv(file_energies)
        merged = energies.merge(estimator,how='inner', left_on='wt', right_on='wt_pdb')
        merged.to_csv('../../output/design/' + designval + '_design.csv', index=None)
 

def mutate (pdb, mutations, pdb_dir, out_dir, repair):

    foldx.create_individual_list(mutations ,'./')
    foldx.run_build_model(pdb, 'individual_list.txt', pdb_dir, out_dir)
    
    foldx.rename_buildmodel_files(pdb[:-4], out_dir, './individual_list.txt')
