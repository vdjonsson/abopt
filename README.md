# abopt. 

AbOpt is a tool for antibody optimization. This computational method jointly constrains the antibody design space and solves for antibody cocktails to address virus escape. Our contributions are: a new statistical model, the saturated LASSO (satlasso) for feature selection and a new combinatorial optimization algorithm that inputs fitness landscapes and solves for antibody cocktails that can address escape through virus mutations. 

## Disclaimer: this code is work in progress currently being cleaned up. ##


# NeurIPS LMRL talk 
https://drive.google.com/file/d/1Zm_ei3fueVl2_HlRLixcX6dNxwPcy2FU/view

# Install 

Dependencies: wget, biopandas, scikit-learn, satlasso 


pip install wget
pip install biopandas
pip install scikit-learn
pip install satlasso 

# Manual

## antibody pipeline 

Typing `antibody_pipeline` will produce the following output: 


    antibody_pipeline v 0.1 
    antibody_pipeline <CMD>

    align         Pseudo align antibody sequences
    distance      Calculate Levenshtein distance between antibody sequences
    estimator     Run estimator on antibody sequences
    map           Map estimator FASTA locations to PDB locations 
    constrain     Constrain estimator features
    scan          Mutational scanning on structure 
    mutate        Mutate a molecular structure 
    repair        Repair a molecular structure
    epitope       Finds molecular structure binding regions
    energy        Run energy minimization functions 
    merge         Merges several batch runs 
    cocktail      Solves for antibody combinations to target virus escape mutants
    version       Prints version information

The usage commands are:

### align

`antibody_pipeline align` pseudo aligns antibody sequences. The arguments for the align command are:

> `input: sequence heavy and light chains `

> `output: matrix of aligned sequences`

### distance
`antibody_pipeline distance` generates a matrix of Levenshtein distances between antibodies. 

    usage: DISTANCE [-h] --filepath FILEP --filename FILEN [--o OUT]

    Calculate pairwise Levenshtein (edit) distance for provided sequences.

    optional arguments:
      -h, --help        show this help message and exit
      --filepath FILEP  Filepath of dataframe
      --filename FILEN  Filename of dataframe (do not include file type extension)
      --o OUT           Output directory for program output

Note: the file contained at filepath/filename.csv must be able to be loaded into a pandas DataFrame with the following required columns:

> `"heavy_chain": column of amino acid sequences for the heavy chain of the antibody`

> `"light_chain": column of amino acid sequences for the light chain of the antibody`

> `"antibody_id": column of identifying string or number for each antibody`

### estimator
`antibody_pipeline estimator` generates an estimator given IC50 data. The arguments for the estimator command are:

    usage: ESTIMATOR [-h] --filepath FILEP --filename FILEN [--s SAT] [--m] --y Y
                 [--cv CV] [--t TRNSFM] --l1 LMBD1 --l2 LMBD2 --l3 LMBD3
                 [--n N_LMBDS] [--r RANGE] [--o OUT]

    Parse amino acid sequence and run SatLasso for variable selection.

    optional arguments:
      -h, --help        show this help message and exit
      --filepath FILEP  Filepath of dataframe
      --filename FILEN  Filename of dataframe (do not include file type extension)
      --s SAT           Saturation value to use for SatLasso(CV): can be float or
                    {"max", "mode"}
      --m               Use flag to map coefficients back to individual amino acid
                    sequences
      --y Y             Name of column for y values
      --cv CV           Cross-validation value: use int for SatLassoCV with
                    specified number of folds; otherwise SatLasso (no CV) used
      --t TRNSFM        Transform to use on y values to feed into cross-
                    validation: can be {"ln", "log10", "norm"}
      --l1 LMBD1        Lambda 1 value; or if using CV, start value for lambda 1
      --l2 LMBD2        Lambda 2 value; or if using CV, start value for lambda 2
      --l3 LMBD3        Lambda 3 value; or if using CV, start value for lambda 3
      --n N_LMBDS       Number of lambdas to use for grid search in CV; ignored if
                    not using CV
      --r RANGE         Range to use for grid search of lambdas in CV; ignored if
                    not using CV
      --o OUT           Output directory for program output

Note: the file contained at filepath/filename.csv must be able to be loaded into a pandas DataFrame with the following required columns:

> `"sequence": column of amino acid sequences to use for one-hot encoding`

The following columns are required only if the -m flag is used:

> `"heavy_chain_aligned": column of aligned amino acid sequences for the heavy chain of the antibody`

> `"light_chain_aligned": column of aligned amino acid sequences for the light chain of the antibody`

> `"antibody_id": column of identifying string or number for each antibody`

### map

`antibody_pipeline map` maps coefficients for a specific antibody from FASTA locations (indexed from 0) to provided PDB locations. The arguments for the map command are:

    usage: MAP [-h] --filepath FILEP --filename FILEN --pdb_filepath FILEP
               --pdb_filename FILEN --name NAME [--o OUT]

    Map locations indexed from 0 to provided PDB locations

    optional arguments:
      -h, --help            show this help message and exit
      --filepath FILEP      Filepath of dataframe
      --filename FILEN      Filename of dataframe (do not include file type
                        extension)
      --pdb_filepath FILEP  Filepath of dataframe mapping FASTA locations to PDB
                        locations
      --pdb_filename FILEN  Filename of dataframe mapping FASTA locations to PDB
                        locations (do not include file type extension)
      --name NAME           Name of antibody in dataframe to map to PDB locations
      --o OUT               Output directory for program output

Note: the file contained at filepath/filename.csv must be able to be loaded into a pandas DataFrame with the following required columns:

> `"antibody_id": column of identifying string or number for each antibody`

> `"location": column of location for each amino acid in the sequence`

> `"chain": column of indicating chain for each amino acid in the sequence`

> `"coefficient": column of coefficients for each amino acid location`

Note: the file contained at pdb_filepath/pdb_filename.csv must be able to be loaded into a pandas DataFrame with the following required columns:

> `"pdb_location": column of PDB location for each amino acid in the sequence`

> `"fasta_location": column of FASTA location for each amino acid in the sequence (must match location in above dataframe)`

> `"chain": column of indicating chain for each amino acid in the sequence`

### constrain
`antibody_pipeline constrain` constrains the features of an estimator. The arguments for the estimator command are:
 constrain(estimator_file, antibody, cutoff , chain, top):

> `input: filename str of path and filename of estimator`

> `input: antibody str of antibody name`

> `input: cutoff array of [cutoffmin, cutoffmax] minimum and maximum coefficient value`

> `input optional: chain str of chain name of the molecular structure to constrain`

> `input optional: top int of number of locations to consider`

> `output: DataFrame of constrained features based on antibody specfic PDB locations`

### mutate
`antibody_pipeline mutate` mutates an antibody given a list of mutations, and generates a structure for each mutation. 

> `input: filename str of molecular structure to mutate`

> `input: repair bool True if structure to mutate requires repair`

> `input: chain str of the chain to repair`

> `input: array list of mutations WT (wild type) amino acid, chain, PDB location, MUT (mutation) amino acid`

> `output: PDB file of molecular structure mutated`

### scan 
`antibody_pipeline scan` mutational scanning on a subset of locations or a chain.

> `input: filename of molecular structures to scan`

> `output: dG of unfolding from molecular structure scanned`

### repair
`antibody_pipeline repair` repairs an antibody using FoldX RepairPDB. 

> `input: array of filenames of molecular structures to repair`

> `input: indir input directory`

> `input: outdir output directory`

> `output: PDB file of molecular structure mutated`


### epitope
`antibody_pipeline epitope` finds epitope, given an amino acid location for a molecular structure.

> `input: molecular structure`

> `input: list of specific locations that includes:  WT (wild type) amino acid, chain, PDB location`

> `input: list of chains to scan for epitopes`

> `output: list of epitopes in the form of chain, PDB locations`

### energy        
`antibody_pipeline energy` generates a matrix of folding energies based on running energy minimization for mutational scanning. The arguments for the energy command are:

> `input: binding ddG or coupling dddG`

> `input: array of molecular structures`

> `input: list of specific locations that includes:  WT (wild type) amino acid, chain, PDB location, MUT (mutation) amino acid`

> `output: matrix of binding energies, or coupling energies`

### merge 
`antibody_pipeline merge` merges energy landscape data for multiple structures

> `input: binding or coupling`

> `input: array of molecular structures`

> `output: merged matrix of binding energies, or coupling energies`

### cocktail
`antibody_pipeline cocktail ` generates combinations of antibodies that are optimal for mutant viruses 

    usage: COCKTAIL[-h] --filepath FILEP --filename FILEN --p P_UNMANAGED --l
                 LMBD [--o OUT]

    Uses convex combinatorial optimization to compute optimal antibody cocktail given fitness lansdacpes and desired mutation coverage.

    optional arguments:
      -h, --help        show this help message and exit
      --filepath FILEP  Filepath of fitness matrix dataframe
      --filename FILEN  Filename of fitness matrix dataframe (do not include file
                    type extension)
      --p P_UNMANAGED   Maximum proportion of viruses not covered by the antibody
                    cocktail
      --l LMBD          Lambda value to use for penalization of number of
                    antibodies chosen
      --o OUT           Output directory for program output
             
### version      
`antibody_pipeline version` prints version number.

### cite 
`antibody_pipeline cite` prints citation information. 
