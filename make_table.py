#!/usr/bin/python

import argparse
import math
import sys
import pandas as pd
import numpy as np
import os
from pathlib import PurePath

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
        if '.pkl' in file:
            df = pd.read_pickle(path)
            df['ChEBI'] = df['ChEBI'].astype(int)
            df = df.set_index('ChEBI')
        else:
            df = pd.read_csv(path, sep='\t', header=None, names=['ChEBI', 'Info'], index_col='ChEBI', dtype={"ChEBI": "int", "Info": "str"})
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

    #gather properties
    data = import_properties()

    for result in results:
        term = PurePath(result).parts[1].split('_ChEBI_IDs.tsv')[0]
        print('making table for %s' % term)

        # import results
        df = pd.read_csv(result, sep = '\t', names=['ChEBI', 'Publication'], dtype={"ChEBI": "int", "Publication": "str"})
        df_results = df.groupby(by=['ChEBI']).count().rename(columns={"Publication": "Count"})

        # make table
        table = make_table(data, df_results)

        # perform normalization
        table.loc[:,"TFIDF"] = table["Count"].astype(float)*table["idf"].astype(float)
        table.loc[:,"TFIDF"] = table.loc[:,"TFIDF"].round(decimals=0).astype(int)
        print(table)

        # write table to file
        write_to_file(table, term)

if __name__ == '__main__':
    main()
