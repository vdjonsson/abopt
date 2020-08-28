import pandas as pd
import seaborn as sb 
import os 
import antibody_pipeline as ap
import energy as e 
    

def mutate(filename, chain, mutations):
    
    # input: filename str of molecular structure to mutate
    # input: repair bool True if structure to mutate requires repair
    # input: chain str of the chain to repair
    # input: array list of mutations WT (wild type) amino acid, chain, PDB location, MUT (mutation) amino acid
    return 


filename = '../../output/estimator/NeutSeqData_VH3-53_66_aligned_mapped_coefficients.csv'

file_ab_locations = '../../data/location/C105_locations.csv'
pdb = '6XCM_Repair'
pdb_less = '6XCM_Repair_less_virus_Repair'
pdbdir = '../../data/pdb/C105/repaired/'
outdir = '../../output/ddg/'


# Constrain estimator based on cutoffs

# ap.constrain(constraintype ='estimator', constrainfile=filename, antibody='C105', cutoff = [-0.4, 0.1], top=10)

# Scan the antibody based on estimator locations 

file_estimator = '../../output/constrain/C105_estimator.csv'
estimator = pd.read_csv(file_estimator).wt_pdb_foldx
# ap.scan (scantype='location', scanvalues = estimator, scanmolecule= 'virus', antibody ='C105', pdbfile = pdb +'.pdb', pdbdir=pdbdir, outdir=outdir)
# ap.scan (scantype='location', scanvalues = estimator, scanmolecule= 'ab', antibody='C105', pdbfile = pdb_less +'.pdb', pdbdir=pdbdir, outdir=outdir)


# Find ddG bind after antibody scanning

indir = '../../data/ddg/C105/'
outdir = '../../output/ddg/'
# ap.energy (antibody='C105', pdb=pdb, pdb_less=pdb_less, scantype = 'ab', energy_type='ddgbind', indir = indir, outdir=outdir)

# Constrain locations where ddG < 0 and some other constraints  for mutation
filename= outdir + 'ddgbind_C105_ab_scanning.txt'
condir = '../../output/constrain/'

ap.constrain (constraintype ='energy', constrainfile=filename, antibody='C105', cutoff = [-1e3, 0.4], top=10)
ap.design (designtype='antibody design', designval = 'C105',file_estimator= condir + 'C105_estimator.csv', file_energies= condir + 'C105_energies.csv')

