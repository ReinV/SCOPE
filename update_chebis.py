#!/usr/bin/python

import networkx
import obonet
import csv
import os.path
import urllib
from urllib.request import urlopen
from bs4 import BeautifulSoup
import datetime
import json
import time
import sys
import pandas as pd


def should_files_be_updated():
    '''
    This function retrieves the timestamp of the latest ChEBI ontology, and the timestamp from a file in in the files folder.
    The timestamps are compared, and if the ontology timestamp is more recent than the file timestamp, a TRUE is returend and the files should be updated.
    '''
    # Constants
    URL = 'http://ftp.ebi.ac.uk/pub/databases/chebi/ontology/'
    ONTOLOGY_COLUMN = 0
    DATE_COLUMN = 1
    DAY_POS = 0
    MONTH_POS = 1
    YEAR_POS = 2
    ONTOLOGY_NAME = 'chebi.obo'
    DATE_SYNTAX = '%d-%b-%Y'
    FILE_TO_CHECK = 'files/ChEBI2Names.tsv'

    # Get ChEBI ontolgoy page into df
    html = urlopen(URL).read()
    soup = BeautifulSoup(html, features="html.parser")
    soup.get_text()
    df = pd.DataFrame([x.split() for x in soup.get_text().split('\n')])

    # Get date from df
    date = df[df[ONTOLOGY_COLUMN] == ONTOLOGY_NAME].values[0][DATE_COLUMN]

    # Create timestamp from date and file in the files folder
    date_object = datetime.datetime.strptime(date, DATE_SYNTAX)
    timestamp_ontology = datetime.datetime.timestamp(date_object)
    timestamp_file = os.path.getmtime(FILE_TO_CHECK)

    return timestamp_ontology > timestamp_file

def get_ontology():
    '''
    This function imports the latest ChEBI ontology.
    '''
    URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
    graph = obonet.read_obo(URL)
    # graph = obonet.read_obo('chebi.obo')
    return graph

def get_mass(node, graph):
    '''
    This function retrieves the mass of a molecule from the ontology.
    '''
    mass = "-"
    try:
        for string in graph.nodes[node]['property_value']:
            if 'mass' in string and 'monoisotopicmass' not in string:
                mass = string.split('\"')[1] # wat te doen met een 0 ?
    except:
        pass
    return mass

def get_smile(node, graph):
    '''
    This function retrieves Smiles from the ontology.
    '''
    smile = ''
    try:
        for value in graph.nodes[node]['property_value']:
            if 'smile' in value:
                smile = value.split('\"')[1]
    except:
        pass
    return smile

def get_relations(nodes, graph):
    '''
    This function recieves a node in the graph, and for that node, it retrieves every parent with 'is_a' relationships.
    It returns a dictionary with those parents.
    '''
    parent_to_key = dict()

    for node in nodes:
        for child, parent, key in graph.out_edges(node, keys=True):
            if key == 'is_a':
                try:
                    parent_to_key[parent]
                except:
                    parent_to_key[parent] = key

    return parent_to_key

def get_superterms(id, graph):
    '''
    This function recieves a ChEBI ID and the latest ChEBI ontology.
    It returns a list of parent ChEBI ids with 'is_a' relationship, up to the root.
    '''
    list_relations = []
    nodes = [id]
    end = False

    while end == False:
        parent_to_key = get_relations(nodes, graph) # get the 'is a' relationships
        if len(parent_to_key) == 0: #if there are no 'is a' relationships, end the search
            end = True
        else:
            nodes = []
            for parent in parent_to_key.keys():
                nodes.append(parent)
                new_id = parent.split(":")[1]
                list_relations.append(new_id)

    return list_relations

def get_new_names(graph, file):
    '''
    This function recieves the latest ChEBI ontolog, and the file containing the ChEBI names.
    It returns a dictionary with names of new molecules not yet stored in the file.
    '''
    id_to_info = read_file(file)
    id_to_name = dict()
    for key, data in graph.nodes(data=True):
        id = key.split(":")[1]
        try:
            id_to_info[id]
        except:
            name = data.get('name')
            id_to_name[id] = name
    return id_to_name

def get_new_smiles(graph, file):
    '''
    This function recieves the latest ChEBI ontology, and the file containing the smiles.
    It returns a dictionary with smiles of new molecules not yet stored in the file.
    '''
    id_to_info= read_file(file)
    id_to_smile = dict()

    for key in graph.nodes():
        id = key.split(":")[1]
        try:
            id_to_info[id]
        except:
            smile = get_smile(key, graph)
            if smile != '':
                id_to_smile[id] = smile
    return id_to_smile

def perform_task(string_of_smiles, SLEEP_TIME, TIMEOUT):
    '''
    '''
    MODEL_ID = 4 # AlogPS3.0 model id
    url = 'https://ochem.eu/modelservice/postModel.do?'
    values = {'modelId': MODEL_ID, 'mol': string_of_smiles}

    # Start task
    connection = False
    while connection == False:
        try:
            data = urllib.parse.urlencode(values)
            data = data.encode('ascii')
            response = urllib.request.urlopen(url, data, timeout=TIMEOUT)
            json_data = json.loads(response.read())
            if json_data["taskId"] == 0:
                connection = False
                print('task failed, trying again in %d seconds' % SLEEP_TIME)
                time.sleep(SLEEP_TIME)
            else:
                connection = True
        except:
            print('connection failed, trying again in %d seconds' % SLEEP_TIME)
            time.sleep(SLEEP_TIME) # in seconds

    return json_data["taskId"]

def download_results(task_id, SLEEP_TIME, TIMEOUT):
    '''
    This function recieves an OCHEM task id, and uses the task id to download the task output.
    '''
    status = 'pending'
    request_url = 'http://ochem.eu/modelservice/fetchModel.do?taskId='+str(task_id)

    while status == 'pending':
        connection = False
        try:
            response = urllib.request.urlopen(request_url, timeout=TIMEOUT)
            task_output = json.loads(response.read())
            id = task_output['taskId']
            connection = True
            status = task_output['status']
        except:
            if connection == True and id == 0:
                status = 'success'
            else:
                print('Task pending, request repeated in %d seconds' % SLEEP_TIME)
                time.sleep(SLEEP_TIME)
        if connection == True:
            if status == 'error':
                if len(task_output['predictions']) == 0:
                    sys.exit('Model predictions failed, script stopped')
            elif status == 'pending':
                print('Model predictions are still pending: request repeated in %d seconds ...' % SLEEP_TIME)
                time.sleep(SLEEP_TIME)


    return task_output

def get_succesful_predictions(output):
    '''
    This function retrieves succesful predictions (predictions that are not 'None').
    '''
    predictions = []
    for pred in output['predictions']:
        if pred != None:
            predictions.append(pred)
    return predictions

def store_predictions(predictions, id_to_logS, id_to_logP, id_to_smile, counter, exceptions):
    '''
    This function retrieves model predictions, connects these to the corresponding ChEBI ID and adds them to the correct dictionary.
    '''
    for output_list in predictions:
        try:
            predictions = output_list['predictions']
            logP = predictions[0]['value']
            logS = predictions[1]['value']
        except:
            exceptions += 1
            logP = '-'
            logS = '-'
        id = list(id_to_smile.keys())[counter]
        id_to_logP[id] = logP
        id_to_logS[id] = logS
        counter += 1
    return id_to_logS, id_to_logP, counter, exceptions

def get_new_predictions(id_to_smile):
    '''
    This function recieves a dictionary with new smiles. It returns two dictionaries with predicted logP and logS values.
    Smiles are put in a string in batches of 200, and added to an url ...
    With this url using the OCHEM API, a task is created where the model AlogPS3.0 is applied ...
    After the task is completed, the output is downloaded using the task_id ...
    Sometimes, output is invalid for a certain SMILE, resulting in a task error. The task needs to be resubmitted for the remaining SMILES.
    '''

    TIMEOUT = 100
    SLEEP_TIME = 10
    BATCH_LENGTH = 200

    counter = 0
    exceptions = 0
    id_to_logP = dict()
    id_to_logS = dict()

    print("Total new predictions: %d" % len(id_to_smile))

    for i in range(0, len(id_to_smile), BATCH_LENGTH):
        print('Getting predictions %d to %d' % (i, i+BATCH_LENGTH))
        difference = BATCH_LENGTH
        end = i + BATCH_LENGTH
        while difference > 0 and counter < len(id_to_smile):
            # get smiles
            string_of_smiles = ''
            smiles = list(id_to_smile.values())[i:end]
            for smile in smiles:
                string_of_smiles += str(smile) + '$$$$'

            # perform_task
            task_id = perform_task(string_of_smiles, SLEEP_TIME, TIMEOUT)
            task_results = download_results(task_id, SLEEP_TIME, TIMEOUT)
            predictions = get_succesful_predictions(task_results)
            id_to_logS, id_to_logP, counter, exceptions = store_predictions(predictions, id_to_logS, id_to_logP, id_to_smile, counter, exceptions)

            difference = difference - len(predictions)
            if difference != 0:
                i = i + len(predictions)

        time.sleep(SLEEP_TIME)
    print('Got prediction errors for %d of %d predictions in total' % (exceptions, counter))

    return id_to_logP, id_to_logS

def read_file(file):
    '''
    This function reads a file and returns a dictionary of the CHEBI ID's.
    If the file does not exists, the file is made and an empty dictionary is returned.
    '''
    id_to_info = dict()

    if os.path.exists(file):
        f = open(file, 'r')
        lines = f.readlines()

        for line in lines:
            line_to_list = line.split('\t')
            id = line_to_list[0]
            info = line_to_list[1].strip() # nodig?
            id_to_info[id] = info
    else:
        f = open(file, 'w') # make file

    return id_to_info

def rewrite_file_to_pkl(graph, file):
    '''
    This function recieves the ontology graph, and the output file path.
    It retrieves hierarchly classes for every ChEBI identifier in the ontology.
    The resulting dataframe is written to a .pkl file.
    '''
    ids_list = []
    info_list = []
    for key in graph.nodes():
        id = key.split(":")[1]
        info = get_superterms(key, graph)
        ids_list.append(id)
        info_list.append(info)
    df = pd.DataFrame({'ChEBI': ids_list, 'Info': info_list})
    df.to_pickle(file)
    print('%s updated' % file)

def rewrite_file(graph, file):
    '''
    This function recieves the latest ChEBI ontology, and the file that will be overwritten.
    Depending on the file argument, it will get the corresponding information (mass, superterms), and writes this to the file.
    '''
    with open(file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for key in graph.nodes():
            id = key.split(":")[1]
            if file == 'files/ChEBI2Mass.tsv':
                info = get_mass(key, graph)
            elif file == 'files/ChEBI2Class.tsv':
                info = get_superterms(key, graph)
            writer.writerow([id, info])
    print('%s updated' % file)

    return

def update_file(id_to_info, file):
    '''
    This function recieves a dictionary and a file.
    Information in the dictionary will be added to that file.
    '''
    with open(file, 'a', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for id in id_to_info:
            info = id_to_info[id]
            writer.writerow([id, info])
    print('%s updated' % file)
    return

def main():

    if should_files_be_updated():
        print('Files need updating ...')

        # Import ontology
        print('Importing ontology')
        graph = get_ontology()

        # Retrieve information from the ontology
        print('Retrieving info from ontology')
        id_to_name = get_new_names(graph, file='files/ChEBI2Names.tsv')
        id_to_smile = get_new_smiles(graph, file='files/ChEBI2Smiles.tsv')
        id_to_logP, id_to_logS = get_new_predictions(id_to_smile)

        # Update the files with new information
        print('Updating files ...')
        update_file(id_to_name, file='files/ChEBI2Names.tsv')
        update_file(id_to_smile, file='files/ChEBI2Smiles.tsv')
        update_file(id_to_logP, file='files/ChEBI2logP.tsv')
        update_file(id_to_logS, file='files/ChEBI2logS.tsv')

        # Some files we rewrite in stead of updating
        rewrite_file(graph, 'files/ChEBI2Mass.tsv')
        rewrite_file_to_pkl(graph, 'files/ChEBI2Class.pkl')

    else:
        print('Files are up-to-date')

if __name__ == '__main__':
    main()
