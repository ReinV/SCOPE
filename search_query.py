#!/usr/bin/python

import urllib.request
import urllib.parse
import json
import time
import csv
import argparse
import datetime


def construct_url(query, pageSize, cursorMark):
    """
    This function constructs a url that searches through the europepmc website and returns the results in json format.
    """

    query_string = urllib.parse.quote(query)

    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + query_string + "&pageSize=" + str(pageSize) + "&resultType=lite&format=JSON&cursorMark=" + str(cursorMark)

    return url


def get_data(url):
    """
    This function recieves a url and returns the json data as a dictionary or list.
    If the connection fails within 5 seconds, it waits 10 minutes to try again, repeated till succes.
    """
    connection = False
    while connection == False:
        try:
            response = urllib.request.urlopen(url, timeout=5)
            data = json.loads(response.read())
            connection = True
        except:
            # if host kills connection, wait 10 minutes and try again (because connection remains False)
            current_time = datetime.datetime.now()
            print('at time %s connection failed, waiting 10 minutes to retry connection...' % current_time.strftime("%H:%M:%S"))
            time.sleep(600)
    return data


def find_publications_with_tmt(data, dict):
    """
    This function looks for publications with text mined terms.
    It returns a dictionary with the ID's of those publications as keys and the source as value.
    """
    for publication in data['resultList']['result']:
        if publication['hasTextMinedTerms'] == 'Y':
            dict[publication['id']] = publication['source']

    return(dict)


def search_publications(query, pageSize):
    """
    This function searches publications with the query and returns a dictionary with publication ID's that contained text mined terms.
    CursorMark is used to go through the search result pages.
    For every page, new ID's plus sources are added to the publications dictionary.
    """
    cursorMark = "*"
    query_data = get_data(construct_url(query, pageSize, cursorMark))
    nextcursormark = query_data['nextCursorMark']
    publications = dict()
    publications = find_publications_with_tmt(query_data, publications)

    while cursorMark != nextcursormark:
        cursorMark = nextcursormark
        query_data = get_data(construct_url(query, pageSize, cursorMark))
        nextcursormark = query_data['nextCursorMark']
        publications = find_publications_with_tmt(query_data, publications)
        time.sleep(5) # be easy on server

    return publications


def get_annotations(publications):
    """
    This function searches through the publications with text mined terms for annotations of type 'Chemicals'.
    From the ChEBI urls, the ChEBI ID's are extracted and returned as values with the publication ID's as keys in a dictionary.
    """

    chebi_dict = dict()
    count = 0
    for pub_id in publications.keys():
        count += 1
        if count%100 == 0:
            print("getting annotations: %d publications of %d total" % (count, len(publications.keys())))
            time.sleep(5) # be easy on server

        source = publications[pub_id]
        url = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=" + str(source) + ":" + str(pub_id) + "&format=JSON"

        annotations_json = get_data(url)
        for element in annotations_json[0]['annotations']:
            if element['type'] == 'Chemicals':
                chebi_url = element['tags'][0]['uri']
                chebi_id = chebi_url.split('_')[1]
                try: #check if publication ID is already in dictionary, if not, then a list as value needs to be created for potentially multiple chebi id's.
                    chebi_dict[pub_id]
                    chebi_dict[pub_id].append(chebi_id)
                except:
                    chebi_dict[pub_id] = [chebi_id]

    return chebi_dict

def write_results(dict, term, query):
    """
    This function writes the ChEBI urls and publication ID's in a seperate csv file.
    Additionally, metadata is written to a text file in the metadata folder
    """
    file = 'results/'+str(term)+'_ChEBI_IDs.tsv'
    count = 0
    uniques = set()

    with open(file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for pub_id in dict.keys():
            for chebi_id in dict[pub_id]:
                count += 1
                uniques.add(chebi_id)
                writer.writerow([chebi_id, pub_id])

    print('%s query results are written to file' % term)

    current_day = datetime.date.today()
    number_of_papers = len(dict.keys())
    number_of_chemicals = count
    number_of_unique_chemicals = len(uniques)

    file = 'metadata/'+str(term)+'.txt'
    f = open(file, 'w')
    f.write('metadata for %s\n' % term
    + 'query: %s\n' % query
    + 'search date: %s\n' % current_day
    + 'number of papers: %d\n' % number_of_papers
    + 'number of chemicals: %d\n' % number_of_chemicals
    + 'number of unique chemicals: %d (note: not all chemicals can be plotted due to missing logP values)' % number_of_unique_chemicals )

def read_input(file):
    '''
    This function reads the input file and returns the query terms in a dictionary.
    The therm used for saving the file is the string before ", " and the query is the string after the ", "
    Every line in the input file should be a new query.
    '''
    f = open(file, 'r')
    input = f.readlines()
    queries = dict()

    for line in input:
        query_list = line.split(',')
        term = query_list[0]
        query = query_list[1].strip()
        queries[term] = query

    return queries

def parser():
    parser = argparse.ArgumentParser(description='A script that parses JSON data from EPC publications and extracts annotations of type Chemicals')
    parser.add_argument('-i', required=True, metavar='input_file', dest='input_file', help='[i] to select input file from the queries folder')
    parser.add_argument('-p', required=False, metavar='pageSize', dest='pageSize', help='[p] to select pageSize')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    input_file = args.input_file
    #pageSize = int(args.pageSize)
    pageSize = 1000

    queries = read_input(input_file)
    for term in queries.keys():
        query = queries[term]
        publications = search_publications(query, pageSize)
        chebi_dict = get_annotations(publications)
        print('%d publications found with search term %s' % (len(chebi_dict.keys()), term) )
        write_results(chebi_dict, term, query)

if __name__ == '__main__':
    main()
