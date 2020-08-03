#!/usr/bin/python

import argparse
import math
import sys
import pandas as pd
import numpy as np
from os import listdir

def read_searches_by_year():
    '''
    This function reads the latest decade of searches_by_year, and returns a dataframe with ChEBI and Publication identifiers.
    '''
    # path = 'searches_by_year/1990-1999_ChEBI_IDs.tsv'
    path = 'searches_by_year/2010-2019_ChEBI_IDs.tsv'
    df = pd.read_csv(path, sep='\t', dtype={"ChEBI": "str", "Publication": "str"})
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
    This function performs the tfidf normalization ( see https://en.wikipedia.org/wiki/Tf%E2%80%93idf )
    '''
    try:
        npubs = id_to_npubs['Count'][id]
    except:
        npubs = 0
    idf = math.log(N/(1+npubs))
    tfidf = idf * count
    # if id == '15377':
    #     print(id, count, N, npubs, idf, tfidf, tfidf )
    return math.floor(tfidf)

def import_properties():
    '''
    This function reads a file and adds the information + id in a dictionary. This is added to another dictionary as value, and the term that describes the
    information (e.g. "Names") as key.
    '''

    files = ['files/ChEBI2Names.tsv', 'files/ChEBI2Mass.tsv', 'files/ChEBI2logP.tsv', 'files/ChEBI2Class.tsv']
    data = dict()
    for file in files:
        key = file.split('2')[1].split('.')[0]
        df = pd.read_csv(file, sep='\t', header=None, names=['ChEBI', 'Info'], index_col='ChEBI')
        id_to_info = df.to_dict()
        data[key] = id_to_info
    return data

def get_info(data, key, id):
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
        # print(data[key])
        column_list = [get_info(data, key, id) for id in df_results.index.to_list()]
        table.loc[:,key] = column_list
    table = table.replace({'Mass': {"-": np.nan}, 'logP': {"-": np.nan}})
    table = table.dropna()
    return table

def write_to_file(table, term):
    '''
    This function writes the table in a .tsv file and names this file after the (shortest) query search term.
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
        files = listdir(input)
        results = [input+'\\'+file for file in files]
    else:
        sys.exit('Error: please give \'file\' or \'folder\' as input type')

    get searches_by_year and statistics
    print('getting searches by year ...')
    df_sby = read_searches_by_year()
    chebi_to_npubs, N = get_statistics(df_sby)

    for result in results:
        print(result)
        term = result.split('\\')[1].split('_ChEBI_IDs.tsv')[0]
        print('making table for %s' % term)

        # import results
        df = pd.read_csv(result, sep = '\t', names=['ChEBI', 'Publication'])
        df_results = df.groupby(by=['ChEBI']).count().rename(columns={"Publication": "Count"})

        # perform normalization
        df_results.loc[:,"TFIDF"] = [normalize(str(id), count, N, chebi_to_npubs) for id, count in zip(df_results.index, df_results.Count)]
        total_tfidf = df_results.TFIDF.sum()
        total_counts = df_results.Count.sum()
        df_results.loc[:,"TFIDF"] = df_results.TFIDF * (total_counts / total_tfidf)
        df_results = df_results.round()
        df_results = df_results.astype({'TFIDF': 'int32'})

        # gather properties
        data = import_properties()

        # make table
        table = make_table(data, df_results)

        # write table to file
        write_to_file(table, term)

if __name__ == '__main__':
    main()
