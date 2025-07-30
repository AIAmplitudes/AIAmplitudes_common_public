# aiamplitudes_common_public
Common public utils library for the "AI for Scattering Amplitudes" project, located at: https://github.com/AIAmplitudes/AIAmplitudes_common_public/

src/
__init__.py- import-level objects and functions. These include the Phi2 and Phi3 symbols, the front and back spaces and relations in various weights and formats, run-polynomials, and their coeffs.
For example use-cases, please explore the Amplitudes_101 jupyter notebook!

download_data.py - contains the scripts to download a (compressed) version of the data files. Call download_all(DATA_DIR_NAME) to download the files.
This is needed to load the symbols, etc., but does not need to be called often! To download the data in the "public" repo rather than the "Phase 1" repo, 
you will need to generate a github access token that has read access from the "data_public" repo and pass it to download_all(). 
Note that some functions involving front/back spaces etc. involve files that have not yet been made public, so this is necessary in order to explore the 
full functionality of the amplitudes package.

file_readers.py - functions to read the files. convert() is a multipurpose function that is called often, while the others are more specialized to different cases.
Most of these are called within functions and should not need to be called by the typical user.

fbspaces.py - Utils for the front and back spaces. This includes coproducts, lookups, etc.

rels_utils.py- Table of homogeneous relations (as defined in the "Bootstrapping a Stress Tensor..." paper), relation eval functions, etc.- this is called often!

commonclasses.py - Common wrappers used for generating and processing relation-instances.

polynomial_utils.py - functions related to analysis of the "d-ray" polynomials mentioned in the "Recurrent Features of Amplitudes..." paper.

notebook/
Amplitudes_101.ipynb- a notebook that displays some simple exploratory data analysis functions of this package.
