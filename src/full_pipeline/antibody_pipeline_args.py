import sys
from align import run_align
from distance import run_distance
from estimator import run_estimator
from map import run_map
from constrain import run_constrain
from scan import run_scan
from mutate import run_mutate
from repair import run_repair
from epitope import run_epitope
from energy import run_energy
from merge import run_merge
from optimize import run_optimize

def print_menu():
    print(
    """
    antibody_pipeline v 0.1
    antibody_pipeline <CMD>

    align         Pseudo align antibody sequences
    distance      Calculate Levenshtein distance between antibody sequences
    estimator     Run estimator on antibody sequences
    constrain     Constrain estimator features
    scan          Mutational scanning on structure
    mutate        Mutate a molecular structure
    repair        Repair a molecular structure
    epitope       Finds molecular structure binding regions
    energy        Run energy minimization functions
    merge         Merges several batch runs
    optimize      Optimizes antibody combinations to cover mutant viruses
    version       Prints version information
    """
    )

def main():
    if len(sys.argv) == 1:
        print_menu()
        exit()
    elif sys.argv[1].lower() == 'align':
        run_align(sys.argv[2:])
    elif sys.argv[1].lower() == 'distance':
        run_distance(sys.argv[2:])
    elif sys.argv[1].lower() == 'estimator':
        run_estimator(sys.argv[2:])
    elif sys.argv[1].lower() == 'map':
        run_map(sys.argv[2:])
    elif sys.argv[1].lower() == 'constrain':
        run_constrain(sys.argv[2:])
    elif sys.argv[1].lower() == 'scan':
        run_scan(sys.argv[2:])
    elif sys.argv[1].lower() == 'mutate':
        run_mutate(sys.argv[2:])
    elif sys.argv[1].lower() == 'repair':
        run_repair(sys.argv[2:])
    elif sys.argv[1].lower() == 'epitope':
        run_epitope(sys.argv[2:])
    elif sys.argv[1].lower() == 'energy':
        run_energy(sys.argv[2:])
    elif sys.argv[1].lower() == 'merge':
        run_merge(sys.argv[2:])
    elif sys.argv[1].lower() == 'optimize':
        run_optimize(sys.argv[2:])
    elif sys.argv[1].lower() == 'version':
        print('Version 1.0.0') # TO DO: print version information
    else:
        print('Incorrect usage. Correct usage shown below: ')
        print_menu()
        exit()
        

if __name__ == "__main__":
    main()
