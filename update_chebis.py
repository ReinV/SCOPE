#!/usr/bin/python

import networkx
import obonet
import csv
import os.path

def return_latest_ontology():
    '''
    This function imports the latest updated version of the ChEBI ontology, and returns the version number and ontology.
    '''
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
    # file = open('files/chebi_180.obo', encoding = 'utf8')
    graph = obonet.read_obo(url)
    # file.close()

    # Mapping from term ID to name
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    version = graph.graph['data-version']
    return version, graph, id_to_name

def return_current_version():
    '''
    This function opens the ChEBI files with id's and names, and returns the version number used to update this file.
    '''
    file = open('files/ontology_version.txt', 'r')
    version = file.read()
    return version

def return_archived_ontology(version):
    '''
    This function returns an archived ontology based on the version number.
    '''
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel' + version + '/ontology/chebi.obo'
    graph = obonet.read_obo(url)
    # id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    return graph

def show_updates(graph_new, graph_old):
    '''
    This function compares two ontologies and returnes the difference in nodes and edges (chemicals and relations).
    '''
    difference_nodes = len(graph_new) - len(graph_old)
    difference_edges = graph_new.number_of_edges() - graph_old.number_of_edges()
    message = 'Newly updated ChEBI ontology contains %d new chemicals and %d new relations' % (difference_nodes, difference_edges)
    return message

def get_mass(node, graph):
    '''
    This function retrieves the mass of a molecule from the ontology.
    '''
    mass = "-"
    try:
        for value in graph.node[node]['property_value']:
            if 'mass' in value and 'monoisotopicmass' not in value:
                mass = value.split('\"')[1] # wat te doen met een 0 ?
    except:
        pass
    return mass

def get_smiles(node, graph):
    '''
    This function retrieves Smiles from the ontology.
    '''
    smile = ''
    try:
        for value in graph.node[node]['property_value']:
            if 'smile' in value:
                smile = value.split('\"')[1]
    except:
        pass
    return smile

def get_relations(nodes, graph):
    parent_to_key = dict()

    for node in nodes:
        for child, parent, key in graph.out_edges(node, keys=True):
            if key == 'is_a':
                try:
                    parent_to_key[parent]
                except:
                    parent_to_key[parent] = key

    return parent_to_key

def get_superterms(id, graph, id_to_name):
    '''
    This function return the is_a terms for a given ChEBI ID.
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
                list_relations.append(id_to_name[parent])

    return list_relations

def update_version_number(number):
    '''
    This function updates the ontology version text file with the version number of the ontology by which the files have been updated.
    '''
    file = open('files/ontology_version.txt', 'w')
    file.write(number)
    return

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

def update_smile(file, graph, id_to_name):
    '''
    This function writes old and new id's with their smile to a .tsv file, and the new smiles will be written to a seperate text file.
    The new id's will be written to the file after the old id's, so that the order is similar to the new smiles text file.
    '''

    id_to_smile = read_file(file)
    new_smiles = dict()
    for key in graph.nodes():
        id = key.split(":")[1]
        try:
            id_to_smile[id]
        except:
            new_smiles[id] = get_smiles(key, graph)

    with open(file, 'w', newline='', encoding="utf-8") as tsvfile: # first write old smiles to smiles file, and the new smiles to the smiles file
        writer = csv.writer(tsvfile, delimiter = '\t')
        for id in id_to_smile.keys():
            smile = id_to_smile[id]
            if smile != '': # make sure no empty smiles are in the file
                writer.writerow([id, smile])
        for id in new_smiles.keys():
            smile = new_smiles[id]
            if smile != '': # make sure no empty smiles are in the file
                writer.writerow([id, smile])
    tsvfile.close()

    f = open('files/new_smiles.txt', 'w') # then add the new smiles to a text file (in the same order)
    for id in new_smiles.keys():
        smile = new_smiles[id]
        if smile != '': # no empty smiles in the file
            f.write(smile+'\n')
    f.close()

def update_file(file, graph, id_to_name):
    '''
    This function recieves the file path, the corresponding file content in a dictionary, and the latest ontology.
    The keys in the latest ontology are CHEBI IDs. Every CHEBI ID is tested in the dictionary to determine if its present in the file.
    If it's not present, the CHEBI ID and its information (smile, name, or superterms) is added to the file.
    '''
    with open(file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for key in graph.nodes():
            id = key.split(":")[1]
            if file == 'files/ChEBI2Names.tsv':
                info = id_to_name[key]
                # info = get_name(key, graph)
            elif file == 'files/ChEBI2Superterms.tsv':
                info = get_superterms(key, graph, id_to_name)
            elif file == 'files/ChEBI2Mass.tsv':
                info = get_mass(key, graph)
            writer.writerow([id, info])

def main():
    files = ['files/ChEBI2Names.tsv','files/ChEBI2Smiles.tsv', 'files/ChEBI2Superterms.tsv', 'files/ChEBI2Mass.tsv']
    current_version = return_current_version()
    latest_version, graph, id_to_name = return_latest_ontology() # graph = ontology

    if current_version == latest_version:
        print('files are up-to-date')
    else:
        print('files need updating')
        graph_old = return_archived_ontology(current_version)
        # updates = show_updates(graph, graph_old)
        # print(updates)

        for file in files:
            # if file == 'files/ChEBI2Smiles.tsv':
            #     update_smile(file, graph, id_to_name)
            if file == 'files/ChEBI2Superterms.tsv':
                update_file(file, graph, id_to_name)
            print('%s updated' % file)

        update_version_number(latest_version)

if __name__ == '__main__':
    main()
