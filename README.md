# AIAA_common
Common utils library.

To install locally, run "pip install -e ." in a clean conda environment.
I am in the process of publishing to pypi!

src/
download_data.py - contains the scripts to download a (compressed) version of lance's files. Call download_all(your_dir) to download the files.
This is useful for EDA and generators, but should not be called every time- for large-scale generation, it may be useful to put these files in a Docker with the module.

file_readers.py - functions to read the files. convert() is a multipurpose function that is called often, while the others are more specialized to different cases.
Most of these are called within functions and should not need to be called by the typical user.

__init__.py- import-level objects. These include the Phi2 and Phi3 symbols, the front and back spaces and relations in various weights and formats, run-polynomials, and their coeffs.
jupyter-notebook coming soon!

fbspaces.py - Utils for the front and back spaces. This includes coproducts, lookups, etc.


