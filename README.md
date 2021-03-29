# abopt. 

AbOpt is a tool for antibody optimization. This computational method jointly constrains the antibody design space and solves for antibody cocktails to address virus escape. Our contributions are: a new statistical model, the saturated LASSO (satlasso) for feature selection and a new combinatorial optimization algorithm that inputs fitness landscapes and solves for antibody cocktails that can address escape through virus mutations. 

## Disclaimer: this code is work in progress currently being cleaned up. ##


# NeurIPS LMRL talk 
https://drive.google.com/file/d/1Zm_ei3fueVl2_HlRLixcX6dNxwPcy2FU/view

# Install 

Dependencies: wget, biopandas, scikit-learn, satlasso 

To install necessary packages: 
>   ` pip install wget `
>   ` pip install biopandas `
>   ` pip install scikit-learn`
>   ` pip install satlasso `

# Manual

## antibody pipeline 

Typing `abopt` will produce the following output: 


    abopt v 0.1 
    abopt <CMD>

    distance      Calculate Levenshtein distance between antibody sequences
    estimator     Run estimator on antibody sequences
    map           Map estimator FASTA locations to PDB locations 
    constrain     Constrain estimator features
    scan          Mutational scanning on structure 
    mutate        Mutate a molecular structure 
    repair        Repair a molecular structure
    epitope       Finds molecular structure binding regions
    energy        Run energy minimization functions 
    merge         Merges fitness landscape data  
    cocktail      Solves for antibody combinations to target virus escape mutants
    version       Prints version information

The usage commands are:

### distance
`abopt distance` generates a matrix of Levenshtein distances between antibodies. 

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
`abopt estimator` generates an estimator given IC50 data. The arguments for the estimator command are:

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

`abopt map` maps coefficients for a specific antibody from FASTA locations (indexed from 0) to provided PDB locations. The arguments for the map command are:

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
`abopt constrain` constrains the features of an estimator. The arguments for the estimator command are:
 constrain(estimator_file, antibody, cutoff , chain, top):

    usage: CONSTRAIN [-h] --filename FILEN --antibody AB [--cutoff CUTOFF CUTOFF] [--chain CH] [--topk K] [--o OUT]

    Constrains the features of an estimator by enforcing cutoff for coefficients, using top - k features, or
    constraining to certain chain.

    optional arguments:
      -h, --help            show this help message and exit
      --filename FILEN      Filepath/filename of estimator
      --antibody AB         Name of antibody
      --cutoff CUTOFF CUTOFF
                            Cutoff array of [cutoffmin, cutoffmax] minimum and maximum coefficient value, please
                            specify in order: MIN MAX
      --chain CH            Chain name of molecular structure to constrain
      --topk K              Top number of locations to consider
      --o OUT               Output directory for program output

### mutate

`abopt mutate` mutates an a structure given a list of mutations, and generates a structure for each mutation/mutation list. 

    usage: MUTATE [-h] --filename FILEN --chain CHAIN [--mlist_filename MFILEN] [--llist_filename LFILEN] [--repair]
                  [--o OUT]

    Mutate a molecular structure.

    optional arguments:
      -h, --help            show this help message and exit
      --filename FILEN      Filepath/filename of molecular structures to mutate
      --chain CHAIN         Chain on structure to mutate used when mutating locations
      --mlist_filename MFILEN
                            Filepath/filename for list of mutations
      --llist_filename LFILEN
                            Filepath/filename for list of locations to mutate
      --repair              Flag to use if structure(s) to mutate require(s) repair after mutation
      --o OUT               Output directory for program output
      

Note: the files at mlist_filename (optional) and llist_filename (optional) must be comma-delimited files of mutations or locations. For example:

> `TH28I, YH58F`

for mlist_filename; OR

> `28-58, 75`

for llist_filename.

### scan 
`abopt scan` mutational scanning on a subset of locations or a chain.

    usage: SCAN [-h] [--backend PKG] --filenames FILENS [FILENS ...] --scantype TYPE [--chain]
                [--llist_filename LFILEN] [--o OUT]

    Performs mutational scanning on a subset of locations or a chain.

    optional arguments:
      -h, --help            show this help message and exit
      --backend PKG         Software program to use as backend for mutational scanning (i.e. FoldX, Rosetta)
      --filenames FILENS [FILENS ...]
                            Filenames of molecular structures to scan
      --scantype TYPE       Type of structure to scan ("antibody" or "virus")
      --chain               Flag to indicate to scan the entire chain
      --llist_filename LFILEN
                            Filename of location ranges of pdb locations to scan, comma delimited in file
      --o OUT               Output directory for program output
      
Note: the file at llist_filename must be a comma-delimited file of locations. For example:

> `400-403,420-423`

### repair
`abopt repair` repairs an antibody for downstream analysis. 

    usage: REPAIR [-h] [--backend PKG] --filenames FILENS [FILENS ...] --toolpath PATH [--o OUT]

    Repairs an antibody using software tool for repairing structures (default: FoldX RepairPDB).

    optional arguments:
        -h, --help            show this help message and exit
        --backend PKG         Software program to use as backend for repair (i.e. FoldX, Rosetta)
        --filenames FILENS [FILENS ...]
                          Filenames of molecular structures to repair
        --toolpath PATH       Pathname of tool used for repair
        --o OUT               Output directory for program output

### epitope
`abopt epitope` finds epitope, given an amino acid location for a molecular structure.

    usage: EPITOPE [-h] --filename FILEN --pdb_filename PDBF --chains CH [CH ...]
                   [--o OUT]

    Calculate epitopes of given structure

    optional arguments:
      -h, --help            show this help message and exit
      --filename FILEN      Filename of specific locations that includes: WT amino
                            acid, chain, PDB location
      --pdb_filename PDBF   Filename including path for molecular structure (PDB
                            file)
      --chains CH [CH ...]  List of chains to scan for epitopes
      --o OUT               Output directory for program output

### energy        
`abopt energy` generates a matrix of folding energies based on running energy minimization for mutational scanning. The arguments for the energy command are:

    usage: ENERGY [-h] --filenames FILENS [FILENS ...] --binding BIND --coupling COUPL [--llist_filename LFILEN]
              [--o OUT]

    Generates a matrix of folding energies based on running energy minimization for mutational scanning.

    optional arguments:
      -h, --help            show this help message and exit
      --filenames FILENS [FILENS ...]
                            Filenames including path of dG unfold files or ddG binding files used for calculation
      --binding BIND        String value to calculate binding ddG
      --coupling COUPL      String value to calculate coupling dddGs
      --llist_filename LFILEN
                            Filename of location ranges of pdb locations to constrain energy calculations, comma
                            delimited in file
      --o OUT               Output directory for program output

Note: the file at llist_filename must be a comma-delimited file of locations. For example:

> `250-300, 400-410`

### merge 
`abopt merge` merges energy fitness landscape data 

    usage: MERGE [-h] --filenames FILENS [FILENS ...] [--norm NORM] [--o OUT]

    Merges energy landscape data for multiple structures.

    optional arguments:
      -h, --help            show this help message and exit
      --filenames FILENS [FILENS ...]
                            Filenames including path of ddG binding for merging
      --norm NORM           Normalization method for binding energies, if any
      --o OUT               Output directory for program output

### cocktail
`abopt cocktail ` generates antibody cocktails that optimize virus mutation coverage and antibody neutralization given fitness landscape data 

    usage: COCKTAIL [-h] --filepath FILEP --filename FILEN --virus_filepath VFILEP --virus_filename VFILEN --p
                    P_UNMANAGED --g1 GAM1 --g2 GAM2 [--o OUT]

    Uses convex combinatorial optimization to compute optimal antibody cocktail for provided virus mutants.

    optional arguments:
      -h, --help            show this help message and exit
      --filepath FILEP      Filepath of fitness matrix dataframe
      --filename FILEN      Filename of fitness matrix dataframe (do not include file type extension)
      --virus_filepath VFILEP
                            Filepath of fitness matrix dataframe for virus mutants
      --virus_filename VFILEN
                            Filename of fitness matrix dataframe (do not include file type extension) for virus
                            mutants
      --p P_UNMANAGED       Maximum proportion of viruses not covered by the antibody cocktail
      --g1 GAM1             Gamma 1 value to use for penalization of number of antibodies chosen
      --g2 GAM2             Gamma 2 value to use for weighting infectivity of virus mutants
      --o OUT               Output directory for program output
             
### version      
`abopt version` prints version number.

### cite 
`abopt cite` prints citation information. 
