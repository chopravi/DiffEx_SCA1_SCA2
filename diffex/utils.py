from itertools import combinations
import pandas as pd
import numpy as np
import scipy.stats as stats
import random

from diffex import constants

IUPHAR_Channels_names = constants.IUPHAR_Channels_names

def clean_dataframe(df):
    
    """ Return a cleaned dataframe with NaN rows removed and duplicate
        fold change measurements averaged """
    
    # Select all rows from the df that don't have NA
    clean_df = df.loc[df['Gene name'].notnull(), :]
    # Select only rows with Gene names that are duplicated
    dup_df = clean_df[clean_df.duplicated(subset='Gene name',keep=False)]
    dup_df = dup_df.sort_values(by=['Gene name'])
    try: # won't work if no duplicates to average
        # Average duplicate fold change measurements
        dup_df = dup_df.groupby('Gene name',as_index=False).mean()
        dup_df = dup_df.round(3)
    except:
        print(f'No duplicated gene names in dataset for df with column 1: {df.columns[1]}')
        pass
    # Drop rows from the original dataframe that are in the duplicate df
    cond = clean_df['Gene name'].isin(dup_df['Gene name'])
    clean_df.drop(clean_df[cond].index, inplace = True)
    clean_df = clean_df.append(dup_df)
    clean_df = clean_df.reset_index(drop=True)
    
    del dup_df
    return clean_df

def get_combinations(names_list, k):
    
    """ Return a list of unique combinations (each element is a tuple)
        from the names_list """
    
    return list(combinations(names_list, k))

def find_pairwise_overlaps(dfs_dict):
    
    """ Return a dataframe with a column holding the overlapping
        set of 'Gene name's for each unique pair of models in the
        dfs_dict """

    model_pairs = get_combinations(names_list=list(dfs_dict.keys()), k=2)
    overlaps_dict = {}
    
    for combi_tuple in model_pairs:
        # create a name to be used for this combination's 
        # dataframe column
        combi_name = '-'.join(combi_tuple)
        # find overlap between the two model gene name columns
        df_1 = dfs_dict[combi_tuple[0]]
        df_2 = dfs_dict[combi_tuple[1]]
        overlap_df = pd.merge(df_1, df_2, on='Gene name')

        overlaps_dict[combi_name] = overlap_df['Gene name']        
        
    overlaps_df = pd.DataFrame(overlaps_dict)
    return overlaps_df
    
def find_triplet_overlaps(dfs_dict):
    
    """ Return a dataframe with a column holding the overlapping
        set of 'Gene name's for each unique unique group of three 
        models in the dfs_dict """

    model_trips = get_combinations(names_list=list(dfs_dict.keys()), k=3)
    overlaps_dict = {}
    
    for combi_tuple in model_trips:
        # create a name to be used for this combination's 
        # dataframe column
        combi_name = '-'.join(combi_tuple)
        # find overlap between the two model gene name columns
        df_1 = dfs_dict[combi_tuple[0]]
        df_2 = dfs_dict[combi_tuple[1]]
        df_3 = dfs_dict[combi_tuple[2]]
        
        overlap_df = pd.merge(df_1, df_2, on='Gene name')
        overlap_df = pd.merge(overlap_df, df_3, on='Gene name')

        overlaps_dict[combi_name] = overlap_df['Gene name']        
        
    overlaps_df = pd.DataFrame(overlaps_dict)
    
    return overlaps_df

def find_quad_overlaps(dfs_dict):
    
    """ Return a dataframe with a column holding the overlapping
        set of 'Gene name's across all models in dfs_dict """

    model_trips = get_combinations(names_list=list(dfs_dict.keys()), k=4)
    overlaps_dict = {}
    
    for combi_tuple in model_trips:
        # create a name to be used for this combination's 
        # dataframe column
        combi_name = '-'.join(combi_tuple)
        # find overlap between the two model gene name columns
        df_1 = dfs_dict[combi_tuple[0]]
        df_2 = dfs_dict[combi_tuple[1]]
        df_3 = dfs_dict[combi_tuple[2]]
        df_4 = dfs_dict[combi_tuple[3]]
        
        overlap_df = pd.merge(df_1, df_2, on='Gene name')
        overlap_df = pd.merge(overlap_df, df_3, on='Gene name')
        overlap_df = pd.merge(overlap_df, df_4, on='Gene name')

        overlaps_dict[combi_name] = overlap_df['Gene name']        
        
    overlaps_df = pd.DataFrame(overlaps_dict)
    
    return overlaps_df

def hypergeometric_test(dfs_dict):
    
    """ Return nothing. Run a hypergeometric test to determine 
        significance of overlaps between dysregulated channels in
        each pair of models in dfs_dict. Print p-values for likelihood
        of channel overlap between each unique pair of models. """
    
    pairwise_overlaps_df = find_pairwise_overlaps(dfs_dict)

    for pair_name in pairwise_overlaps_df.columns:
        # Define total number of genes and channels overlapping between current two models
        overlapping_genes = pairwise_overlaps_df.loc[:, pair_name].dropna()
        overlapping_channels = pd.Series(list(set(overlapping_genes).intersection(set(IUPHAR_Channels_names))))
        total_channel_overlaps = len(overlapping_channels)

        # Define the total number of channels in the IUPHAR channel database
        IUPHAR_chan_num = len(IUPHAR_Channels_names)

        # Find the names of the two models under consideration
        overlap_model_names = pair_name.split('-')

        # Find the total number of channels dysregulated in each model's dataset
        model_1_genes = dfs_dict[overlap_model_names[0]].loc[:, 'Gene name']
        model_1_channels = pd.Series(list(set(model_1_genes).intersection(set(IUPHAR_Channels_names))))
        total_model_1_channels = len(model_1_channels)

        model_2_genes = dfs_dict[overlap_model_names[1]].loc[:, 'Gene name']
        model_2_channels = pd.Series(list(set(model_2_genes).intersection(set(IUPHAR_Channels_names))))
        total_model_2_channels = len(model_2_channels)

        pairwise_overlap_p_value = 1-stats.hypergeom.cdf(total_channel_overlaps,
                                                         IUPHAR_chan_num,
                                                         total_model_1_channels,
                                                         total_model_2_channels)

        print(f'{overlap_model_names[0]} and {overlap_model_names[1]} p-value={pairwise_overlap_p_value}')

def get_n_channels_dict(dfs_dict):
    
    """ Return a dict holding the number of 
        channels dysregulated in each model """
    
    n_channels_dict = {}
    
    for model_name, model_df in dfs_dict.items():
        
        gene_names = model_df.loc[:, 'Gene name']
        channel_names = pd.Series(list(set(gene_names).intersection(set(IUPHAR_Channels_names))))
        n_channels = len(channel_names)
        
        n_channels_dict[model_name] = n_channels
        
    return n_channels_dict

def set_channels_df(dfs_dict, filename):
    # Save and get a csv with each model as column and its respective
    # list of dysregulated channels along rows
    
    model_channels_df_dict = {}
    for model_name, model_df in dfs_dict.items():
        
        gene_names = model_df.loc[:, 'Gene name']
        channel_names = pd.Series(list(set(gene_names).intersection(set(IUPHAR_Channels_names))))
        model_channels_df_dict[model_name] = channel_names
    
    model_channels_df = pd.DataFrame(model_channels_df_dict)
    if filename:
        model_channels_df.to_csv(filename)
    else:
        pass
    
    return model_channels_df, model_channels_df_dict

def simulate_channel_dfs(dfs_dict, n_channels_dict):
    
    """ Return a dict of dataframes simulated according to the data
        in dfs_dict. For each model, choose a rand set of n channels
        from the IUHPAR database where n = number of channels dysregulated
        in that model. n is defined in n_channels_dict """
    
    # for each model, choose a random set of of channels from IUPHAR database
    # of size however meany channels are in that models original data    
    sim_dfs_dict = {}

    for model_name, model_df in dfs_dict.items():
        # get the number of channels dysregulated in this model
        n_channels = n_channels_dict[model_name]
        # Sample the IUPHAR channels and make a simulated dataframe
        # for this model's dysregulated Genes
        sim_channel_names = IUPHAR_Channels_names.sample(n_channels)
        sim_model_df = pd.DataFrame({'Gene name': sim_channel_names})

        # Add this simulated model df to the sim_dfs_dict
        sim_dfs_dict[model_name] = sim_model_df
        
    return sim_dfs_dict

def simulate_triplet_overlaps(dfs_dict, n_runs, observed_overlap_n, n_channels_dict):
    
    """ Return a tuple: (triplet_overlaps, n_succeses) where triplet_overlaps is a
        list of total overlaps found for each simulation of length n_runs and
        n_successes is the number of runs with >= observed_overlap_n simulated
        overlaps. """
    
    triplet_overlaps = []
    n_successes = 0
    
    for i in range(0, n_runs):
        
        # Simulate a dictionary of dataframes, one dataframe
        # for each model
        sim_channel_dfs_dict = simulate_channel_dfs(dfs_dict, n_channels_dict)
        sim_triplet_overlaps_df = find_triplet_overlaps(sim_channel_dfs_dict)
        triplet_overlaps_list = [sim_triplet_overlaps_df.loc[:, column].dropna() for column in sim_triplet_overlaps_df.columns]
        triplet_overlaps_names = pd.concat(triplet_overlaps_list).unique() # returns a series object
        
        triplet_overlaps_num = len(triplet_overlaps_names)
        print(f'Run number {i+1} of {n_runs}. Found {triplet_overlaps_num} triple overlapping channels',
              end='\r')
        
        triplet_overlaps.append(triplet_overlaps_num)
        
        if triplet_overlaps_num >= observed_overlap_n:
            n_successes += 1
        else:
            pass
        
    return triplet_overlaps, n_successes

def simulate_quad_overlaps(dfs_dict, n_runs, observed_overlap_n, n_channels_dict):
    
    """ Return a tuple: (quad_overlaps, n_succeses) where quad_overlaps is a
        list of total overlaps found for each simulation of length n_runs and
        n_successes is the number of runs with >= observed_overlap_n simulated
        overlaps. """
    
    quad_overlaps = []
    n_successes = 0
    
    for i in range(0, n_runs):
        
        # Simulate a dictionary of dataframes, one dataframe
        # for each model
        sim_channel_dfs_dict = simulate_channel_dfs(dfs_dict, n_channels_dict)
        sim_quad_overlaps_df = find_quad_overlaps(sim_channel_dfs_dict)
        quad_overlaps_list = [sim_quad_overlaps_df.loc[:, column].dropna() for column in sim_quad_overlaps_df.columns]
        quad_overlaps_names = pd.concat(quad_overlaps_list).unique() # returns a series object
        
        quad_overlaps_num = len(quad_overlaps_names)
        print(f'Run number {i+1} of {n_runs}. Found {quad_overlaps_num} quad overlapping channels',
              end='\r')
        
        quad_overlaps.append(quad_overlaps_num)
        
        if quad_overlaps_num >= observed_overlap_n:
            n_successes += 1
        else:
            pass
        
    return quad_overlaps, n_successes

def drop_non_channels(overlaps_df, filename):    
    
    """ Return the overlap dataframe with all channels dropped
        and index reset. Save the df as a csv with the filename
        passed this function. """
    
    df = overlaps_df
    
    channels_df_dict = {}
    for column in df.columns:
        # For each set of overlaps, drop all the gene names that are not
        # channels. They are replaced by NaNs.
        channels_bool = df.loc[:, column].isin(IUPHAR_Channels_names)
        channels_df_dict[column] = df.loc[channels_bool, column]

    channels_df = pd.DataFrame(channels_df_dict)
    
    clean_channels_df = channels_df.reset_index(drop=True).copy()    
    for column in channels_df.columns:
        # Set all of the rows in this column to NaN so they can be replaced
        # by lists of channel names in each overlap.
        clean_channels_df.loc[:, column] = np.NaN
        channel_names = list(channels_df.loc[:, column].dropna())
        # Put the list of channels in the overlap's row. Save the df
        clean_channels_df.loc[0:len(channel_names)-1, column] = channel_names
        clean_channels_df.to_csv(filename)
        
    return clean_channels_df