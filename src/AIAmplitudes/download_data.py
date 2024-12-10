from __future__ import annotations

import tempfile
import tarfile
from pathlib import Path

import requests

#__all__ = ("data_path", "download_all", "known_files")

DIR = Path(__file__).parent.resolve()

zipurl = "https://github.com/AIAmplitudes/data_public"

def _cache_path(cache_dir: str | None = None) -> Path:
    if cache_dir is None:
        ampdir = Path.home() / ".local" / "AIAmplitudesData"
        ampdir.mkdir(exist_ok=True, parents=True)
        return ampdir
    return Path(cache_dir)

local_default = _cache_path(None)

def download_all(cache_dir: str | None = None) -> None:
    local_dir = _cache_path(cache_dir)
    with tempfile.TemporaryFile() as f:
        with requests.get(zipurl, stream=True) as r:
            r.raise_for_status()
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

        tar = tarfile.open("test.tar")
        tar.getmembers()
        with tarfile.TarFile(f) as z:
            print(z.namelist())
            #for n in z.namelist():
            #    if "src/AIAmplitudes/data/" in n:
            #        z.extract(n, str(local_dir / str(n.split("/")[-1])))
