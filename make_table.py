#!/usr/bin/python


import argparse
import math
import sys
import pandas as pd
import numpy as np
import os
from pathlib import PurePath
# from os import listdir

def read_searches_by_year():
    '''
    This function reads the summarized searches by decade files in the searches_by_year folder,
    and returns a dataframe with ChEBI and Publication identifiers.
    '''
    folder = 'searches_by_year/'
    files = [file for file in os.listdir(folder) if '.tsv' in file]
    df = pd.DataFrame(columns=["ChEBI", "Publication"])
    for file in files:
        path = folder+file
        df_new = pd.read_csv(path, sep='\t', dtype={"ChEBI": "str", "Publication": "str"})
        df = pd.concat([df, df_new])
    return df

def get_statistics(df):
    '''
    - Extract total amount of unique publiacations.
    - Construct dictionary with ChEBI ID to the amount of unique publications containing that ChEBI ID.
    '''
    N = len(set(df.Publication))
    df = df.groupby(by=['ChEBI']).count()
    df = df.rename(columns={"Publication": "Count"})
    chebi_to_count = df.to_dict()
    return chebi_to_count, N

def normalize(id, count, N, id_to_npubs):
    '''
    This function performs the tfidf normalization based on the number of publications a ChEBI identifers can be found in
    ( see https://en.wikipedia.org/wiki/Tf%E2%80%93idf ).
    '''
    try:
        npubs = id_to_npubs['Count'][id]
    except:
        npubs = 0
    idf = math.log(N/(1+npubs))
    tfidf = idf * count
    return math.floor(tfidf)

def import_properties():
    '''
    This function reads a file and adds the information + id in a dictionary. This is added to another dictionary as value, and the term that describes the
    information (e.g. "Names") as key.
    '''
    FOLDER = 'files'
    files = os.listdir(FOLDER)
    data = dict()
    for file in files:
        path = os.path.join(FOLDER, file)
        key = file.split('2')[1].split('_')[0]
        if '.pkl' in file :
            df = pd.read_pickle(path)
            df['ChEBI'] = df['ChEBI'].astype(int)
            df = df.set_index('ChEBI')

        else:
            df = pd.read_csv(path, sep='\t', header=None, names=['ChEBI', 'Info'], index_col='ChEBI')
        id_to_info = df.to_dict()
        data[key] = id_to_info
    return data

def get_info(data, key, id):
    '''
    This function recieves:
        - a ChEBI identifier
        - a dictionary key indicating the type of information
        - data dictionary containing information of every ChEBI identifier
    The function returns specific information retrieved from the data dictionary (e.g. Name or Mass)
    '''
    try:
        info = data[key]['Info'][id]
    except:
        info = float('NaN')
    return info

def make_table(data, df_results):
    '''
    This function recieves the data of all the ChEBI files in the files folder and the ids of the query search.
    It returns a dictionary of the query ids and their properties from the ChEBI files, if those properties are there.
    e.g. if there is no logP value, the id is not added to the dictionary that is returned.
    '''
    table = df_results
    # create column lists
    for key in data.keys():
        column_list = [get_info(data, key, id) for id in df_results.index.to_list()]

        table.loc[:,key] = column_list

    table = table.replace({'Mass': {"-": np.nan}, 'logP': {"-": np.nan}})
    table = table.dropna()
    table = table.sort_values(by='Count', ascending=False)
    return table

def write_to_file(table, term):
    '''
    This function writes the table in a .pkl file for easy importation into the visualization script,
    and a .tsv for own inspection. The file is named with the (shortest) query search term.
    '''
    path = 'tables/'+term
    table.to_csv(path+'_table.tsv', sep='\t')
    table.to_pickle(path+'_table.pkl')

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input', dest='input', help='[i] to select input folder or input file from the results folder ')
    parser.add_argument('-t', required=True, metavar='type', dest='type', help='[t] to select type of input: file or folder')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    input_type = args.type
    input = args.input
    if input_type == 'file':
        results = [input]
    elif input_type == 'folder':
        files = os.listdir(input)
        results = [input+'/'+file for file in files]
    else:
        sys.exit('Error: please give \'file\' or \'folder\' as input type')

    # # get searches_by_year and statistics
    print('getting searches by year ...')
    df_sby = read_searches_by_year()
    chebi_to_npubs, N = get_statistics(df_sby)

    # gather properties
    data = import_properties()

    for result in results:
        term = PurePath(result).parts[1].split('_ChEBI_IDs.tsv')[0]
        print('making table for %s' % term)

        # import results
        df = pd.read_csv(result, sep = '\t', names=['ChEBI', 'Publication'], dtype={"ChEBI": "int", "Publication": "str"})
        df_results = df.groupby(by=['ChEBI']).count().rename(columns={"Publication": "Count"})

        # perform normalization
        df_results.loc[:,"TFIDF"] = [normalize(str(id), count, N, chebi_to_npubs) for id, count in zip(df_results.index, df_results.Count)]
        total_tfidf = df_results.TFIDF.sum()
        total_counts = df_results.Count.sum()
        df_results.loc[:,"TFIDF"] = df_results.TFIDF * (total_counts / total_tfidf)
        df_results = df_results.round()
        df_results = df_results.astype({'TFIDF': 'int32'})

        # make table
        table = make_table(data, df_results)

        # write table to file
        write_to_file(table, term)

if __name__ == '__main__':
    main()
