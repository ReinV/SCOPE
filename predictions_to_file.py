#!/usr/bin/python

import argparse
import csv

def read_smiles(file):
    '''
    This function recieves the ChEBI2Smiles file, it reads the file and returns two dictionaries.
    One dictionary has the id's as keys and smiles as values, and the other has positions of the id's with smiles as keys and the id's as values.
    '''
    f = open(file, 'r')
    data = f.readlines()
    position_to_id = dict()

    position = 0
    for line in data:
        id = line.split()[0]
        position_to_id[position] = id
        position += 1

    f.close()
    return position_to_id

def get_predictions(file, position_to_id):
    '''
    This function recieves the model predictions file and the position_to_id dictionary that helps connecting the predictions with the corresponding id.

    Because only newly addes smiles are used for predictions, and these new smiles were added at the bottum of the ChEBI2Smiles.tsv file,
    the while loop starts at the end (last positions) and continues the loop until all the predictions are linked to the correct ChEBI ID.

    Then, the function returns a dictionary with the id's as keys and their predictions as values.
    '''
    f = open(file, 'r')
    predictions = f.readlines()
    id_to_predictions = dict()

    total_length = len(position_to_id.keys())
    new_smiles_length = len(predictions) - 1 # do not count header
    starting_position = total_length - new_smiles_length

    j = 1 # skip header
    for i in range(starting_position,total_length):
        id = position_to_id[i]
        prediction = predictions[j].split(",")
        logP = prediction[0].strip('\"')
        logS = prediction[1].strip()
        logS = logS.strip('\"')
        id_to_predictions[id] = {"logP": logP, "logS": logS}
        j = j + 1

    return id_to_predictions

def write_to_file(id_to_predictions, name):
    '''
    This function recieves the id's and predictions in a dictionary and writes them to the output file.
    '''
    file = 'files/ChEBI2'+str(name)+'.tsv'

    with open(file, 'a+', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for id in id_to_predictions.keys():
            info = id_to_predictions[id][name]
            writer.writerow([id, info])
    return

def parser():
    parser = argparse.ArgumentParser(description='A script that reads model predictions and writes them with corresponding chebi IDs in a .tsv file')
    parser.add_argument('-m', required=False, metavar='modelpredictions', dest='modelpredictions', help='[c] to select modelpredictions file')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    modelpredictions = args.modelpredictions
    chebi2smiles = 'files/ChEBI2Smiles.tsv'

    position_to_id = read_smiles(chebi2smiles)
    id_to_predictions = get_predictions(modelpredictions, position_to_id)

    write_to_file(id_to_predictions, 'logP')
    write_to_file(id_to_predictions, 'logS')

if __name__ == '__main__':
    main()
