import numpy as np
import pandas as pd

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

def one_hot_encode(aa):
    """
    One-hot encode an amino acid according to amino acid list (aalist) specified above.
    
    Parameters
    ----------
    aa : str
        One-character representation of an amino acid
    
    Returns
    ----------
    encoding : ndarray of shape (len_aalist,)
        One-hot encoding of amino acid
    """
    
    if aa not in aalist:
        return [0]*len(aalist)
    else:
        encoding = [0]*len(aalist)
        encoding[aalist.index(aa)] = 1
        return encoding

def seqparser(df, seq_col):
    """
    Parse amino acid sequences in dataframe; create a one-hot encoded matrix of sequences.
    
    Note: amino acid sequences must have the same length in order to be parsed.
    
    Parameters
    ----------
    df : pandas DataFrame
        Dataframe containing amino acid sequences to parse
    seq_col : str, default = 'sequence'
        Column in dataframe with amino acid sequences
    
    Returns
    ----------
    aamatrix : numpy array of shape (n_sequences, (len_sequence * len_aalist))
        One-hot encoded matrix of amino acid sequences in dataframe
    """
    
    aamatrix = np.empty((0, len(df[seq_col][0])*len(aalist)), int)
    for seq in df[seq_col]:
        row = []
        for aa in seq:
            row = row + one_hot_encode(aa)
        aamatrix = np.vstack((aamatrix,row))
    return aamatrix

def create_coef_matrix(coefs):
    """
    Utility function for map_coefs function.
    
    Creates pandas DataFrame with indices : amino acids, columns : positions in sequence /
        and values : coefficients
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV)
    
    Parameters
    ----------
    coefs : pandas DataFrame
        Coefficients with indices : amino acid positions of form <amino acid><position>
    
    Returns
    ----------
    aa_posmap : pandas DataFrame
        Coefficients with indices : amino acids and columns : positions in sequence
    """
    
    aa_posmap = pd.DataFrame(index=aalist)
    for i in range(0, int(len(coefs)/len(aalist))):
        index = coefs.index[i*len(aalist)][1:]
        aa_posmap[index] = coefs.iloc[:,0][i*len(aalist):(i+1)*len(aalist)].values
    return aa_posmap

def map_coefs(df, coefs, heavy_chain_name, light_chain_name, id_col, seq_col):
    """
    Maps aligned sequences of amino acids back to original non-aligned sequences
    
    Note: Intended to be used with coefficient output from regression package /
        (ex. SatLasso, SatLassoCV)
    
    Parameters
    ----------
    df : pandas DataFrame
        Metadata including amino acid sequences for heavy and light chains, identifying name /
            of each heavy/light chain pair
    coefs : pandas DataFrame
        Coefficients with indices : amino acid positions of form <amino acid><position>
    heavy_chain_name : str
        Name of column for heavy chain sequences in df.
    light_chain_name : str
        Name of column for light chain sequences in df.
    id_col : str
        Name of column for identifier for each heavy/light chain sequences /
            (ex. name of antibody)
    
    Returns
    ----------
    coefs_df : pandas DataFrame
        Coefficients mapped back to amino acid in each heavy/light chain sequence with /
            MultiIndex : identifying_name, location, chain, amino_acid and columns: /
            wild_type and coefficient
    """
    
    if heavy_chain_name is not None and light_chain_name is not None:
        coefs_df = pd.DataFrame(columns = [id_col, 'location', 'chain', 'aa', 'wild_type', 'coefficient'])
        coefs_df.set_index([id_col, 'location', 'chain', 'aa'], inplace=True)

        aa_posmap = create_coef_matrix(coefs)

        len_heavy_chain = len(df[heavy_chain_name][0])
        for antibody in df[id_col]:
            sequence = df[heavy_chain_name][df[id_col] == antibody].item()
            pos = 0
            for i in range(0, len(sequence)):
                if sequence[i] in aalist:
                    wt = [False] * len(aalist)
                    wt[aalist.index(sequence[i])] = True
                    coefs_df = coefs_df.append(pd.DataFrame.from_dict({id_col: [antibody]*len(aalist), 'location': [str(pos)]*len(aalist), 'chain': ['H']*len(aalist), 'aa': aalist, 'wild_type': wt, 'coefficient': aa_posmap[str(i)]}, orient = 'columns').set_index([id_col, 'location', 'chain', 'aa']))
                    pos = pos+1

        for antibody in df[id_col]:
            sequence = df[light_chain_name][df[id_col] == antibody].item()
            pos = 0
            for i in range(0, len(sequence)):
                if sequence[i] in aalist:
                    wt = [False] * len(aalist)
                    wt[aalist.index(sequence[i])] = True
                    coefs_df = coefs_df.append(pd.DataFrame.from_dict({id_col: [antibody]*len(aalist), 'location': [str(pos)]*len(aalist), 'chain': ['L']*len(aalist), 'aa': aalist, 'wild_type': wt, 'coefficient': aa_posmap[str(i+len_heavy_chain)]}, orient = 'columns').set_index([id_col, 'location', 'chain', 'aa']))
                    pos = pos+1
    
    else:
        coefs_df = pd.DataFrame(columns = [id_col, 'location', 'aa', 'wild_type', 'coefficient'])
        coefs_df.set_index([id_col, 'location', 'aa'], inplace=True)

        aa_posmap = create_coef_matrix(coefs)
        
        for antibody in df[id_col]:
            sequence = df[seq_col][df[id_col] == antibody].item()
            pos = 0
            for i in range(0, len(sequence)):
                if sequence[i] in aalist:
                    wt = [False] * len(aalist)
                    wt[aalist.index(sequence[i])] = True
                    coefs_df = coefs_df.append(pd.DataFrame.from_dict({id_col: [antibody]*len(aalist), 'location': [str(pos)]*len(aalist), 'aa': aalist, 'wild_type': wt, 'coefficient': aa_posmap[str(i)]}, orient = 'columns').set_index([id_col, 'location', 'aa']))
                    pos = pos+1
                
    return coefs_df
