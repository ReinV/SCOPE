#!/usr/bin/python

import os
import requests
import json 
import shutil
import tqdm

def download_files(file_to_link, folder):
    '''
    This function downloads the missing files from the OSF project (https://osf.io/pvwu2/).
    '''
    for file, link in file_to_link.items():
        r = requests.get(link, stream=True)
        file_size = int(r.headers['Content-Length'])
        desc = 'Downloading %s' % file
        path = os.path.join(folder, file)

        with tqdm.tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc) as r_raw:
            with open(path, "wb") as f:
                shutil.copyfileobj(r_raw, f)

    return

def does_file_exist(file, folder):
    '''
    This function returns boolean if file exist or not.
    If the folder doesn't exist, it will be created.
    '''
    if not os.path.isdir(folder):
        os.mkdir(folder)

    path = os.path.join(folder, file)

    return os.path.exists(path)

def get_response(url):
    '''
    This function returns url response.
    '''
    response = requests.get(url)
    response = json.loads(response.text)
    return response

def get_file_to_link(url, folder):
    '''
    This function returns a dictionary with download links for files that are missing from the repository.
    '''
    response = get_response(url)
    file_to_link = {}
    for object in response["data"]:
        file_name = object["attributes"]["name"]

        if not does_file_exist(file_name, folder):
            download_link = object["links"]["download"]
            file_to_link[file_name] = download_link

    return file_to_link

def main():
    # Define constants
    FILES_FOLDER = 'files'
    SEARCHES_FOLDER = 'searches_by_year'
    STORAGE_ID_FILES = '611252ba847d1304ca38b4d4'
    STORAGE_ID_SEARCHES = '5f3a6bb93bf2520111a1f53f'
    base_url = 'https://api.osf.io/v2/nodes/pvwu2/files/osfstorage/'

    # Get a list of files from OSF, download the files that are not present in the repository 
    url = base_url + STORAGE_ID_FILES
    file_to_link = get_file_to_link(url, FILES_FOLDER)
    download_files(file_to_link, FILES_FOLDER)

    # Again for searches folder
    url = base_url + STORAGE_ID_SEARCHES
    file_to_link = get_file_to_link(url, SEARCHES_FOLDER)
    download_files(file_to_link, SEARCHES_FOLDER)

if __name__ == '__main__':
    main()
