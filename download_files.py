#!/usr/bin/python

import os
import requests
import shutil
import tqdm

def return_file_info():
    '''
    '''
    file_to_id = {
    'ChEBI2logP': 'zjrmy',
    'ChEBI2logS': '5nuew',
    'ChEBI2Mass': '3s7bn',
    'ChEBI2Class': 'dh7kr',
    'ChEBI2Names': 'gmq9h',
    'ChEBI2Smiles': 'jusvb'
    }
    return file_to_id
# def return_file_info():
#     '''
#     '''
#     file_to_id = {
#     'ChEBI2logP.tsv': 'zjrmy',
#     'ChEBI2logS.tsv': '5nuew',
#     'ChEBI2Mass.tsv': '3s7bn',
#     'ChEBI2Class.pkl': 'dh7kr',
#     'ChEBI2Names.tsv': 'gmq9h',
#     'ChEBI2Smiles.tsv': 'jusvb'
#     }
    # return file_to_id

def return_decade_info():
    '''
    '''
    decade_to_id = {
    '1940-1949_ChEBI_IDS.tsv': 'dz2ks',
    '1950-1959_ChEBI_IDS.tsv': 'rhyad',
    '1960-1969_ChEBI_IDS.tsv': 'v5gsp',
    '1970-1979_ChEBI_IDS.tsv': '9emna',
    '1980-1989_ChEBI_IDS.tsv': 't4ksw',
    '1990-1999_ChEBI_IDS.tsv': 'gy2fw',
    '2000-2009_ChEBI_IDS.tsv': 'rea8g',
    '2010-2019_ChEBI_IDs.tsv': 'cfjde',
    '2020-2029_ChEBI_IDs.tsv': 'qptm4',
    'pre1945_ChEBI_IDs.tsv': '5cnua'
    }

    return decade_to_id

def check_for_files(dict, folder):

    new_dict = {}

    print('Checking for folder: %s' % folder)
    if not os.path.isdir(folder):
        os.mkdir(folder)

    else:
        files = os.listdir(folder)
        for key, id in dict.items():
            if not is_in_files(key, files):
                new_dict[key] = id

    return new_dict

def is_in_files(key, files):
    is_in_files = False
    for file in files:
        if key in file:
            is_in_files = True
    return is_in_files

def check_if_files_exist(dict, folder):
    '''
    This function first checks if the 'Files' and 'Searches by year' folder exist.
    If not, it creates the folders.

    It also deletes entries from the dictionary that are present in the projfect.

    It returns a dictionary retaining only missing files for downloading.
    '''

    print('Checking for folder: %s' % folder)
    if not os.path.isdir(folder):
        os.mkdir(folder)

    new_dict = {}
    for file, id in dict.items():

        if not os.path.exists(folder+file):

            new_dict[file] = id

    return new_dict

def get_filename(r):
    '''
    '''
    info = r.headers['Content-Disposition']
    info_listed = info.split(";")
    file_name = info_listed[1].split("\"")[1]
    return file_name

def download_missing_files(dict, folder, base_url):
    '''
    This function downloads the missing files from the OSF project (https://osf.io/pvwu2/).
    '''

    for file, id in dict.items():
        url = base_url + id + '/'
        r = requests.get(url, stream=True)


        file_name = get_filename(r)
        file_size = int(r.headers['Content-Length'])
        desc = 'Downloading %s' % file_name
        path = folder+file_name

        with tqdm.tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc) as r_raw:
            with open(path, "wb") as f:
                shutil.copyfileobj(r_raw, f)

    return

def main():
    # Define constants
    FILES_FOLDER = 'files/'
    SEARCHES_FOLDER = 'searches_by_year/'
    URL = 'https://osf.io/download/'

    # Get file path + OSF id, hardcoded in this script
    file_to_id = return_file_info()
    decade_to_id = return_decade_info()

    # Go file by file to see if it alread exists in project folder
    # Download missing file from the OSF storage
    file_to_id = check_for_files(file_to_id, FILES_FOLDER)
    download_missing_files(file_to_id, FILES_FOLDER, URL)

    # Again for searches folder
    decade_to_id = check_for_files(decade_to_id, SEARCHES_FOLDER)
    download_missing_files(decade_to_id, SEARCHES_FOLDER, URL)

if __name__ == '__main__':
    main()
