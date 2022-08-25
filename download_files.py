#!/usr/bin/python

import os
import requests
import json 
import shutil
import tqdm
import collections

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


def get_response(url):
    '''
    This function returns url response.
    '''
    response = requests.get(url)
    response = json.loads(response.text)
    return response

def get_osf_to_rel(url):
    '''
    This function ....
    '''
    response = get_response(url)
    ofs_to_rel = collections.defaultdict(dict)

    for object in response["data"]:
        file_name = object["attributes"]["name"]
        download_link = object["links"]["download"]

        if "_" in file_name:
            rel = file_name.split("_")[1].split(".")[0]
            file = file_name.split("_")[0]
        else:
            rel = ""
            file = file_name.split(".")[0]
            
        ofs_to_rel[file]["rel"] = rel
        ofs_to_rel[file]["link"] = download_link
        ofs_to_rel[file]["name"] = file_name

    return ofs_to_rel

def get_repo_to_rel(folder):
    '''
    This function ....
    '''
    files = os.listdir(folder)
    repo_to_rel = collections.defaultdict(dict)

    for file_name in files:
        rel = file_name.split("_")[1].split(".")[0]
        file = file_name.split("_")[0]
        path = os.path.join(folder, file_name)

        repo_to_rel[file]["rel"] = rel
        repo_to_rel[file]["path"] = path
        
    return repo_to_rel

def get_files_to_download(osf_to_rel, repo_to_rel):
    
    files_to_download = {}

    for file in osf_to_rel.keys():
        rel = osf_to_rel[file]['rel']
        link = osf_to_rel[file]['link']
        name = osf_to_rel[file]['name']

        if file in repo_to_rel.keys():

            if rel != repo_to_rel[file]["rel"]:

                files_to_download[name] = link
                path = repo_to_rel[file]["path"]
                os.remove(path)
        else:

            files_to_download[name] = link


    return files_to_download

def main():
    # Define constants
    folder = 'files'
    url = 'https://api.osf.io/v2/nodes/pvwu2/files/osfstorage/611252ba847d1304ca38b4d4/'

    # create if not exists 
    if not os.path.isdir(folder):
        os.mkdir(folder)

    # Get a list of files from OSF, download the files that are not present in the repository or have newer release
    osf_to_rel = get_osf_to_rel(url)
    repo_to_rel = get_repo_to_rel(folder)

    files_to_download = get_files_to_download(osf_to_rel, repo_to_rel)

    download_files(files_to_download, folder)

if __name__ == '__main__':
    main()
