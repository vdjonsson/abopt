import pandas as pd 
import argparse

'''
### merge 
`abopt merge` merges energy landscape data for multiple structures

> `input: files array filenames including pathnames, of ddG binding for merging`
> `input: normalization str normalization method for binding energies, if any`
> `output: merged_raw file with a merged matrix of binding energies, or coupling energies`
> `output: merged_norm file with a merged and normalized matrix of binding energies, or coupling energies`
'''


def setup_parser():
    parser = argparse.ArgumentParser(prog = 'MERGE', description = 'Merges energy landscape data for multiple structures.')
    parser.add_argument('--filenames', nargs = '+', type = str, dest = 'filenames', required = True, help = 'Filenames including path of ddG binding for merging', metavar = 'FILENS')
    parser.add_argument('--norm', type = str, default=None, dest = 'normalization', help = 'Normalization method for binding energies, if any', metavar = 'NORM')
    parser.add_argument('--o', type = str, required = False, default = '', dest = 'output_dir', help = 'Output directory for program output', metavar = 'OUT')
    return parser


def create_output_structure(output_dir):
    if os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir+'output/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/merge/')
        except FileExistsError:
            pass
        try:
            os.mkdir(output_dir+'output/figs/')
        except FileExistsError:
            pass
    else:
        raise KeyError('Output directory '+output_dir+' not found.')


def merge_ab_bloom(ab_name, mut, normalize=True):
    bloom, bloomval = get_bloom_data(normalize=normalize)
    ddg, ddgval = get_ab_ddg_data(ab_name, mut, normalize=normalize)
    merged = bloom.merge(ddg, on='mutation')
    #merged.to_csv(out_tmp + mut+'tmp.csv')                                                                                                
    return merged,bloomval,ddgval

def combine_data_sets(ddg_array):

    merged = pd.DataFrame()
    merged = ddg_array[0]

    for ddg_data in ddg_array[1:]: 
        merged = merged.merge(ddg_data, on='mut') 
    
    merged = merged.drop_duplicates(subset=['mut'], keep='first')
    merged = merged.set_index('mut')
    return merged 


def get_melted_ddgs(ddgs):
    ddgs = ddgs.reset_index()
    ddgs['location'] = ddgs.mut.str[3:-1]
    ddgsm = pd.melt(ddgs, id_vars=['mut', 'location'], value_vars=ddgs.columns[1:-1], var_name = 'ab',value_name='ddg')
    return ddgsm

'''API'''
def combine_antibody_binding_energy_data(abnames, scantype='virus'):
    ddgv = pd.DataFrame()
    for ab in abnames:
        print(ab)
        abddg = get_ab_ddg_data(ab, mutation = None, scantype=scantype, normalize = False)
        ddgv = pd.concat([abddg, ddgv])

    u.write(ddgv, filename='ddgs_virus_scanning')
    return ddgv

def check_args(namespace):
    # TODO
    pass

def run_merge(args):
    parser = setup_parser()
    namespace = parser.parse_args(args)
    check_args(namespace)
    
    create_output_structure(namespace.output_dir)
    
    ###### IMPLEMENTATION #######
       
    # `input: filenames including path of ddG binding for merging, comma delimited` ===> namespace.filenames (this is an array of strings which are each filename)
    # `input: normalization str normalization method for binding energies, if any` ===> namespace.normalization
    
      # output directory: namespace.output_dir (the output should go in namespace.output_dir/output/energy/) and figs in namespace.output_dir/output/figs/
