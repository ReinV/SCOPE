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

def get_latest_release():
    '''
    '''
    URL = 'http://ftp.ebi.ac.uk/pub/databases/chebi/archive/'

    # Get ChEBI ontolgoy page into df
    html = urlopen(URL).read()
    soup = BeautifulSoup(html, features="html.parser")
    soup.get_text()
    df = pd.DataFrame([x.split() for x in soup.get_text().split('\n')])

    release = df.sort_values(by=0, ascending=False).iloc[0,0]
    release = release.replace("/", "")

    return release

def should_files_be_updated(files, latest_release):
    '''
    '''
    update_these_files = []
    for file in files:
        release = file.split(".")[0].split("_")[1]

        if latest_release > release:
            update_these_files.append(file)

    return update_these_files


def get_ontology():
    '''
    This function imports the latest ChEBI ontology.
    '''
    URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
    graph = obonet.read_obo(URL)
    # graph = obonet.read_obo('chebi.obo')
    return graph

def get_mass(graph, file):
    '''
    This function retrieves the mass of a molecule from the ontology.
    '''
    id_to_info = read_file(file)
    id_to_mass = {}
    for node in graph.nodes():
        id = node.split(":")[1]
        try:
            id_to_info[id]
        except:
            mass = "-"
            try:
                for string in graph.nodes[node]['property_value']:
                    if 'mass' in string and 'monoisotopicmass' not in string:
                        mass = string.split('\"')[1] # wat te doen met een 0 ?
            except:
                pass

            id_to_mass[id] = mass
    return id_to_mass

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
    print(url)
    # Start task
    connection = False
    while connection == False:
        try:
            data = urllib.parse.urlencode(values)
            data = data.encode('ascii')
            response = urllib.request.urlopen(url, data, timeout=TIMEOUT)
            json_data = json.loads(response.read())
            print(json_data)
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

    print(request_url)
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
                print(task_output)
                # if len(task_output['predictions']) == 0:
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
    BATCH_LENGTH = 100

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

def rewrite_file_to_pkl(file, graph):
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
            info = get_superterms(key, graph)
            writer.writerow([id, info])
    print('%s updated' % file)

    return

def write_to_file(file, id_to_info):
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

def update_file(file, graph, id_to_smiles, id_to_logP, id_to_logS):
    '''
    '''
    if "Mass" in file:
        id_to_mass = get_mass(graph, file)
        write_to_file(file, id_to_mass)

    if "Names" in file:
        id_to_mass = get_new_names(graph, file)
        write_to_file(file, id_to_mass)

    if "Smiles" in file:
        write_to_file(file, id_to_smiles)

    if "logP" in file:
        write_to_file(file, id_to_logP)

    if "logS" in file:
        write_to_file(file, id_to_logS)

    if "Class" in file:
        rewrite_file_to_pkl(file, graph)

def rename_file(file, latest_release):
    '''
    '''
    p1 = file.split("_")[0]
    p2 = file.split(".")[1]
    new_name = p1+"_"+latest_release+"."+p2

    os.rename(file, new_name)

# def get_latest_release():


def main():
    # Get file names from the files folder
    files = [os.path.join('files', file) for file in os.listdir('files')]

    # Get latest ChEBI ontology release
    latest_release = get_latest_release()

    # Check if files should be updated
    update_these_files = should_files_be_updated(files, latest_release)

    if len(update_these_files) > 0:
        print('Files need updating')

        # Import ontology
        print('Importing ontology...')
        graph = get_ontology()

        print("Updating Files...")

        smiles_file_listed = [file for file in update_these_files if "Smiles" in file]
        logp_file_listed = [file for file in update_these_files if "logP" in file]
        logs_file_listed = [file for file in update_these_files if "logS" in file]

        if len(smiles_file_listed) > 0:
            smiles_file = smiles_file_listed[0]
            id_to_smiles = get_new_smiles(graph, smiles_file)

            if len(logp_file_listed) > 0 or len(logs_file_listed) > 0:
                id_to_logP, id_to_logS = get_new_predictions(id_to_smiles)
            else:
                id_to_logP = None
                id_to_logS = None
        else:
            id_to_logP = None
            id_to_logS = None
            id_to_smiles = None

        for file in update_these_files:
            update_file(file, graph, id_to_smiles, id_to_logP, id_to_logS)
            rename_file(file, latest_release)

    else:
        print('Files are up-to-date')
    #


if __name__ == '__main__':
    main()
