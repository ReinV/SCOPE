#!/usr/bin/python

import os.path
import urllib
from urllib.request import urlopen

def check_and_download():
    '''
    This function first checks if the 'Files' and 'Searches by year' folder exist.
    If not, it creates the folders and downloads the files from the OSF project (https://osf.io/pvwu2/).
    '''
    # Constants
    FILES_FOLDER = 'files/'
    LIST_OF_FILE_ID = ['dh7kr', 'zjrmy', '3s7bn', '5nuew', '5nuew', 'gmq9h', 'jusvb']
    PICKLE_FILE_ID = 'dh7kr'
    PICKLE_FILE_NAME = 'ChEBI2Class.pkl'
    SEARCHES_FOLDER = 'searches_by_year/'
    LIST_OF_SEARCHES_ID = ['qptm4', 'cfjde', 'rea8g']
    URL = 'https://osf.io/download/'

    if not os.path.isdir(FILES_FOLDER):
        os.mkdir(FILES_FOLDER)
        for id in LIST_OF_FILE_ID:
            download_url = URL+id+'/'

            with urllib.request.urlopen(download_url) as response, open(FILES_FOLDER+response.info().get_filename(), 'wb') as out_file:
                print(response.info().get_filename())
                data = response.read()
                out_file.write(data)

    if not os.path.isdir(SEARCHES_FOLDER):
        'Downloading searches by year files, this may take a while to download'
        os.mkdir(SEARCHES_FOLDER)
        for id in LIST_OF_SEARCHES_ID:
            download_url = URL+id+'/'
            with urllib.request.urlopen(download_url) as response, open(SEARCHES_FOLDER+response.info().get_filename(), 'wb') as out_file:
                print(response.info().get_filename())
                data = response.read()
                out_file.write(data)

    return

def main():
    check_and_download()

if __name__ == '__main__':
    main()
