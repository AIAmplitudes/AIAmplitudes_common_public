from __future__ import annotations
import tempfile
import tarfile
from bs4 import BeautifulSoup
import json
import os
import requests
from pathlib import Path

username= None
password = None

################### Download tarballs from git ###############################
public_repo =  "AIAmplitudes/data_public"
private_repo =  "AIAmplitudes/data_private"

def _cache_path(cache_dir: str | None = None) -> Path:
    if cache_dir is None:
        ampdir = Path.home() / ".local" / "AIAmplitudesData"
        try:
            ampdir.mkdir(exist_ok=True, parents=True)
        except:
            print('Could not make dir. Defaulting to local cache')
            os.mkdir('./cache')
            ampdir='./cache'
        return Path(ampdir)
    return Path(cache_dir)

relpath = _cache_path(None)

def get_gitfilenames(the_zipurl):
    soup = BeautifulSoup(requests.get(the_zipurl).text)
    files=[]
    for elem in soup.find_all('script', type='application/json'):
        if ".tar" in elem.text:
            files += [i["name"] for i in json.loads(elem.contents[0])["props"]["initialPayload"]["tree"]["items"]]
    return files

def download_unpack(myfile: str, local_dir: Path):
    with tempfile.TemporaryFile() as f:
        with requests.get(myfile, stream=True) as r:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        f.seek(0)
        with tarfile.open(fileobj=f) as tarf:
            for n in tarf.getnames():
                assert os.path.abspath(os.path.join(local_dir, n)).startswith(str(local_dir))
            tarf.extractall(path=local_dir)
    return

def download_all_public(repo: str = public_repo, cache_dir: str | None = None) -> None:
    #use BS4 to get all the datasets in the public repo
    local_dir = _cache_path(cache_dir)
    if not len(os.listdir(local_dir))==0:
        print("Local cache not empty! Terminating")
        return
    print(f"downloading files from {repo}, unpacking in {local_dir}")
    url=f"https://github.com/{repo}"
    for file in get_gitfilenames(url):
        if not ".tar" in file: continue
        myfile = f"https://raw.githubusercontent.com/{repo}/main/{file}"
        print(f"extracting {myfile}")
        download_unpack(myfile,local_dir)

    #dump all files into the root directory
    for subdir, dirs, files in os.walk(local_dir):
        for file in files:
            os.rename(str(os.path.join(subdir, file)),str(os.path.join(local_dir, file)))

    #delete all subdirs that do not have a file in them
    deleted=set()
    for thisdir, subdirs, files in os.walk(local_dir, topdown=False):
        still_has_subdirs = False
        for subdir in subdirs:
            if os.path.join(thisdir, subdir) not in deleted:
                still_has_subdirs = True
                break
        if not any(files) and not still_has_subdirs:
            os.rmdir(thisdir)
            deleted.add(thisdir)
    return

def download_all_private(mytoken, repo: str = public_repo, cache_dir: str | None = None) -> None:
    #use the github CLI and token access to get datasets in the private repo
    g = Github(token)
    repo_dir = relpath
    git_url = "https://github.com/AIAmplitudes/data_public.git"
    Repo.clone_from(git_url, repo_dir)
    return

def download_all(mytoken = None):
    if mytoken: download_all_private(mytoken)
    download_all_public()
    return
#######################################################################################
