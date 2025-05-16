# script to generate json files of a user-defined set of linear relations matched with a given symbol
# script to evaluate various linear relations satisfied by the symbols of the 3-point form factor of phi2


import numpy as np
import time
import datetime
import re
import itertools
from itertools import permutations
import random
import json
import copy

alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
quad_prefix = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
steinmanns={'a':'d', 'b':'e', 'c':'f', 'd':'aef', 'e':'bdf', 'f':'cde'}

#How many rels of each type to generate? (phase-1 relgen approach)
rels_to_generate= {"first": [[500]*3,[0]*3],
                            "double": [[500]*3,[0]*3],
                            "triple": [[500],[1]],
                            "dihedral": [[500],[1]],
                            "final": [[500]*19+[500]*10,[0]*19+[1]*10],
                            "integral": [[500]*3,[1]*3]
                              }

rels_to_generate_compact_default = {'first': [[500]*3, [0]*3],
                                    'double': [[500]*3, [0]*3],
                                    'triple': [[500], [1]],
                                    'integral': [[500]*3, [1]*3]}


# Removing duplicates
def dropdups(dicts):
    return [dict(t) for t in {frozenset(d.items()) for d in dicts}]

def dropmdups(dicts):
    # also drop scalar multiples
    def normalize(dict):
        #sort keys into alphabetical order
        firstval=sorted(dict.items())[0][1]
        return {k:(int(elem) if (elem := v/firstval).is_integer() else elem) for k,v in dict.items()}
    return dropdups([normalize(dict) for dict in dicts])

##############################################################################################
# HOMOGENOUS LINEAR RELATIONS LOOK-UP TABLES#
##############################################################################################

######### Nonlocal rels #################
# triple-adjacency relation: plus dihedral symmetry; any slot
# integrability relations: any slot
dihedral_table = [list(permutations(alphabet[:3]))[i]+list(permutations(alphabet[3:]))[i]
                  for i in range(len(alphabet))]
cycle_table = [dihedral_table[i] for i in [0, 3, 4]]
flip_table = [dihedral_table[i] for i in [0, 1, 2, 5]]
triple_table = [{'aab': 1, 'abb': 1, 'acb': 1}]
pair_table = [{'ab': 1, 'ac': 1, 'ba': -1, 'ca': -1},  # eq 3.6
                              {'ca': 1, 'cb': 1, 'ac': -1, 'bc': -1},  # eq 3.7
                              {'db': 1, 'dc': -1, 'bd': -1, 'cd': 1, 'ec': 1, 'ea': -1, 'ce': -1,
                               'ae': 1, 'fa': 1, 'fb': -1, 'af': -1, 'bf': 1, 'cb': 2,'bc': -2},
                                 {'ad':1},{'da':1},{'df':1}]

######### Localized rels #################
# first entry condition
first_entry_rel_table = [{'d': 1}, {'e': 1}, {'f': 1}]  # Sec 3.1 (iv)

# double-adjacency condition: plus dihedral symmetry; any slot
double_adjacency_rel_table = [{'de': 1}, {'ad': 1}, {'da': 1}]  # eq 2.19, 2.20

# triple-adjacency relation: plus dihedral symmetry; any slot
triple_adjacency_rel_table = [{'aab': 1, 'abb': 1, 'acb': 1}]  # eq 2.21

# integrability relations: any slot
integral_rel_table = [{'ab': 1, 'ac': 1, 'ba': -1, 'ca': -1},  # eq 3.6
                      {'ca': 1, 'cb': 1, 'ac': -1, 'bc': -1},  # eq 3.7
                      {'db': 1, 'dc': -1, 'bd': -1, 'cd': 1, 'ec': 1, 'ea': -1, 'ce': -1,
                       'ae': 1, 'fa': 1, 'fb': -1, 'af': -1, 'bf': 1, 'cb': 2,
                       'bc': -2}]  # eq 3.8 coeff (8)! Takes longest time.

# multi-final-entries relations: plus dihedral symmetry
# new order: one-term relations, short relations (<=4 terms), long relations (>4 terms).
final_entries_rel_table = [{'a': 1}, {'b': 1}, {'c': 1},  # eq 4.6 (idx: 0-2)
                           {'ad': 1}, {'ed': 1},  # eq 4.7 (1) (idx: 3-4)
                           {'add': 1}, {'abd': 1}, {'ace': 1}, {'ebd': 1}, {'edd': 1},  # eq 4.9 (idx: 5-9)
                           {'addd': 1}, {'abbd': 1}, {'adbd': 1}, {'cbbd': 1},  # eq 4.10 (idx: 10-13)
                           {'ebbd': 1}, {'ebdd': 1}, {'edbd': 1}, {'eddd': 1}, {'fdbd': 1},  # eq 4.11 (idx: 14-18)

                           {'bf': 1, 'bd': -1},  # eq 4.7 (2) (idx: 19)
                           {'cdd': 1, 'cee': 1},  # eq 4.8 (1) (idx: 20)
                           {'ddbd': 1, 'dbdd': -1},  # eq 4.12 (idx: 21)
                           {'cbdd': 1, 'cdbd': -1},  # eq 4.15(1) (idx: 22)
                           {'fbd': 1, 'dbd': -1, 'bdd': 1},  # eq 4.8 (2) (idx: 23)
                           {'bddd': 1, 'faff': 1, 'dbdd': -1, 'eaff': -1, 'fbdd': 1, 'aeee': -1},  # eq 4.13 (idx: 24)
                           {'abdd': 1, 'cddd': -1 / 2, 'dcee': -1 / 2, 'aeee': 1 / 2, 'eaff': 1 / 2, 'faff': -1 / 2,
                            'ecee': 1 / 2},  # eq 4.14 (idx: 25)
                           {'cbdd': 1, 'bfff': -1 / 2, 'dcee': 1 / 2, 'ecee': -1 / 2, 'cddd': 1 / 2, 'dbdd': 1 / 2,
                            'fbdd': -1 / 2},  # eq 4.15(2) (idx: 26)
                           {'cdbd': 1, 'bfff': -1 / 2, 'dcee': 1 / 2, 'ecee': -1 / 2, 'cddd': 1 / 2, 'dbdd': 1 / 2,
                            'fbdd': -1 / 2},  # eq 4.15(3) (idx: 27)
                           {'fbbd': 1, 'dbbd': -1, 'bbdd': 1, 'faff': -1 / 2, 'dbdd': 1 / 2, 'fbdd': -1 / 2,
                            'eaff': 1 / 2, 'aeee': 1 / 2, 'bfff': -1 / 2}]  # eq 4.16 (idx: 28)

# multi-initial-entries relations: plus dihedral symmetry
# new order: one-term relations, short relations (<=4 terms), long relations (>4 terms).
initial_entries_rel_table = [{'ad': 1},
                             {'aad': 1}, {'bcf': 1}, {'bde': 1}, {'bdf': 1}, {'bda': 1}, {'abd': 1},
                             {'cb': 1, 'bc': -1},
                             {'cd': 1, 'bd': -1},
                             {'aaf': 1, 'bbf': 1, 'abf': -1},
                             {'aab': 1, 'aac': 1, 'cca': 1, 'bba': -1, 'aba': -1},
                             {'bba': 1, 'bbc': 1, 'ccb': 1, 'aab': -1, 'abb': -1},
                             {'abc': 1, 'aac': 1, 'bbc': 1, 'cca': 1, 'ccb': 1},
                             {'aac': 1, 'cca': 1, 'bbc': -1, 'ccb': -1, 'afa': 1 / 2, 'aaf': -1 / 2, 'bbf': 1 / 2,
                              'afb': -1 / 2}]

all_rel_table={'first': first_entry_rel_table,
               'double': double_adjacency_rel_table,
               'triple': triple_adjacency_rel_table,
               'integral': integral_rel_table,
               'final': final_entries_rel_table,
               'initial': initial_entries_rel_table}


def trivial_zero_rel_table(format="full"):
    myrel_table = first_entry_rel_table
    slots = [0] * len(first_entry_rel_table)

    if format == "full":
        myrel_table += final_entries_rel_table[:3]
        slots += [-1] * len(final_entries_rel_table[:3])

    steinmanns = get_rel_table_dihedral(double_adjacency_rel_table)

    myrel_table += steinmanns
    slots += [None] * len(steinmanns)
    return myrel_table, slots


def get_image(word, rownum):
    return ''.join([dihedral_table[rownum][idx] for idx in [alphabet.index(l) for l in [*word]]])

def table_image(table):
    return [{get_image(k,ind):v for k,v in rel.items()} for ind in range(len(alphabet)) for rel in table]

def table_to_rels(table):
    tr=table_image(table)
    return [i for n, i in enumerate(tr) if not
               any(set(sorted(i.keys())) == set(sorted(k.keys())) for k in tr[n + 1:])]

pair_rels=table_to_rels(pair_table)
triple_rels=table_to_rels(triple_table)
############################################################################################################
def sumstring(i, mystring, k, v):
    if v == 1:
        if i == 0:
            mystring += f'{k}'
        else:
            mystring += f'+{k}'
    elif v == -1:
        mystring += f'-{k}'
    elif v > 0:
        if i == 0:
            mystring += f'{v}*{k}'
        else:
            mystring += f'+{v}*{k}'
    elif v < 0:
        mystring += f'{v}*{k}'
    else:
        mystring = mystring

    return mystring

def sumdict(k, d1, d2):
    out = {}
    for k in d1 | d2:
        if (k in d1) and (k in d2):
            val = d1[k] + d2[k]
        elif (k in d1):
            val = d1[k]
        elif (k in d2):
            val = d2[k]
        else:
            val = 0
    out[k] = val
    return out


def find_all(a_str, sub):
    for match in re.finditer(re.escape(sub), a_str):
        yield match


    #start = 0
    #while True:
    #    start = a_str.find(sub, start)
    #    if start == -1: return
    #    yield start
    #    start += 1




###################################################################################
########################################################################################################################
# Auxiliary Functions
#######################################################################################################################

#needed
def read_rel_info(rels_to_generate, make_zero_rels=False):
    '''
    Parse a rel_info dict.
    ---------
    INPUTS:
    rels_to_generate: dict; rel type and ID.

    OUTPUTS:
    rels, slots, to_gens, overlaps, relnames: lists; consisting of: relation dicts, slots for rel, num to generate per rel, how much to overlap with sym, name of rel.
    '''

    rels, slots, to_gens, overlaps, relnames = [], [], [], [], []
    for rel_key, rel_info in rels_to_generate.items():
        if rel_key == 'first':
            myrel_table = get_rel_table_dihedral(first_entry_rel_table)
            myslot = 0
        elif rel_key == 'initial':
            myrel_table = get_rel_table_dihedral(initial_entries_rel_table)
            myslot = 0
        elif rel_key == 'double':
            myrel_table = get_rel_table_dihedral(double_adjacency_rel_table)
            myslot = None
        elif rel_key == 'triple':
            myrel_table = get_rel_table_dihedral(triple_adjacency_rel_table)
            myslot = None
        elif rel_key == 'integral':
            myrel_table = get_rel_table_dihedral(integral_rel_table)
            myslot = None
        elif rel_key == 'final':
            myrel_table = get_rel_table_dihedral(final_entries_rel_table)
            myslot = -1
        elif rel_key == 'dihedral':
            myrel_table = [None] * len(rel_info[0])
            myslot = None
        else:
            print("unknown relation!")
            raise ValueError

        for i in range(len(rel_info[0])):
            if (not make_zero_rels) and (rel_info[1][i] == 0): continue
            rels.append(myrel_table[i])
            slots.append(myslot)
            to_gens.append(rel_info[0][i])
            if len(myrel_table[i]) == 1: overlaps.append(0)
            else: overlaps.append(rel_info[1][i])
            relnames.append(f'{rel_key}_{i}')

    return rels, slots, to_gens, overlaps, relnames

def read_allrel_info(rels_to_generate, make_zero_rels=False):
    '''
    Prepare to generate all rels of the types in the rel_info dict- ignore ordering
    ---------
    INPUTS:
    rels_to_generate that does not specify subtype. rels_to_gen is
    {"initial":[1000],
    "final":[2000], etc.}"

    OUTPUTS:
    rels, slots, to_gens, overlaps, relnames: lists; consisting of: relation dicts, slots for rel, num to generate per rel, how much to overlap with sym, name of rel.
    '''

    rels, slots, to_gens, overlaps, relnames = [], [], [], [], []
    for rel_key, rel_info in rels_to_generate.items():
        if rel_key == 'first':
            myrel_table = dropmdups(get_rel_table_dihedral(first_entry_rel_table))
            myslot = 0
        elif rel_key == 'initial':
            myrel_table = dropmdups(get_rel_table_dihedral(initial_entries_rel_table))
            myslot = 0
        elif rel_key == 'double':
            myrel_table = dropmdups(get_rel_table_dihedral(double_adjacency_rel_table))
            myslot = None
        elif rel_key == 'triple':
            myrel_table = dropmdups(get_rel_table_dihedral(triple_adjacency_rel_table))
            myslot = None
        elif rel_key == 'integral':
            myrel_table = dropmdups(get_rel_table_dihedral(integral_rel_table))
            myslot = None
        elif rel_key == 'final':
            myrel_table = dropmdups(get_rel_table_dihedral(final_entries_rel_table))
            myslot = -1
        elif rel_key == 'dihedral':
            myrel_table = [None] * len(rel_info[0])
            myslot = None
        else:
            print("unknown relation!")
            raise ValueError

        #establish a canonical order for the rels of a given type: sort by num letters,
        #then shortest to longest, then alphabetical by keys as list
        myrel_table.sort(key= lambda d:(len(next(iter(d))), (len(d)), list(d.keys())))
        i, mylen=0,0
        for rel in myrel_table:
            if rel_key == "final" and len(next(iter(rel))) != mylen:
                i = 0
                mylen=len(next(iter(rel)))
            if (not make_zero_rels) and (len(rel) == 1): continue

            rels.append(rel)
            slots.append(myslot)
            to_gens.append(rel_info[0][0])
            overlaps.append(0 if (len(rel) == 1) else 1)
            if rel_key == "final":
                relnames.append(f'{rel_key}{mylen}_{i}')
            else:
                relnames.append(f'{rel_key}_{i}')
            i += 1

    return rels, slots, to_gens, overlaps, relnames


def get_coeff_from_word(word, symb):
    '''
    Get the coeff of a given word in a symbol.
    ---------
    INPUTS:
    word: str; a string of letters such as 'aaae'.
    symb: dict; a dictionary with word as key, and coeff as value, e.g., {'aaae':16, 'aaaf':16}.

    OUTPUTS:
    coeff: float; if word does not exist in symb, then return 0.
    '''
    if word in symb:
        return symb[word]
    return 0
def get_word_from_coeff(coeff, symb):
    '''
    Get the words corresponding to a given coeff in a symbol.
    ---------
    INPUTS:
    coeff: float; e.g., 16.
    symb: dict; a dictionary with word as key, and coeff as value, e.g., {'aaae':16, 'aaaf':16}.

    OUTPUTS:
    word: set; a set of words (strings) with the given coeff;
          if coeff does not exist in symb, then return an empty string.
    '''
    if coeff in symb.values():
        word = {i for i in symb if symb[i] == coeff}  # set
        return word
    return str()
def get_dihedral_images(word):
    '''
    Get all the dihedral images of a given word.
    ---------
    INPUTS:
    word: str.

    OUTPUTS:
    dihedral_images: list; each item in the list is a word (str); always has six items.
    '''
    word_idx = [alphabet.index(l) for l in [*word]]
    dihedral_images = [''.join([dihedral_table[row][idx] for idx in word_idx]) for row in range(len(alphabet))]
    return dihedral_images
def get_valid_dihedral_images(word, pruned_symb, badsymb):
    '''
    Get all the dihedral images of a given word that are in a symb and not in a badsymb.
    ---------
    INPUTS:
    word: str.

    OUTPUTS:
    dihedral_images: list; each item in the list is a word (str); always has six items.
    '''
    word_idx = [alphabet.index(l) for l in [*word]]
    dihedral_images = {row: image for row in range(len(alphabet)) if (
        image := ''.join([dihedral_table[row][idx] for idx in word_idx])) in pruned_symb and image not in badsymb}
    return dihedral_images
def get_dihedral_pair(key, goodkeys, symb, type="cycle"):
    '''
    Given a key, the
    ---------
    INPUTS:
    key: str; key in dict.
    key_images. list; the dihedral images of the key.
    goodkeys. set; the allowed keys.
    symb. dict; the symbol to lookup coeffs.
    OUTPUTS:
    pair; dict. A two-term instance of type "cycle" or "flip".
    '''
    # print(goodkeys)
    if type == "cycle":
        indices = set([3, 4])
    elif type == "flip":
        indices = set([1, 2, 5])
    else:
        raise ValueError
    good_inds = list(indices.intersection(set(goodkeys.keys())))
    if len(good_inds) == 0: return None
    # get random key: ind from goodkeys
    ind = random.choice(good_inds)
    pair = {key: [symb[key], 1], goodkeys[ind]: [symb[goodkeys[ind]], -1]}
    return pair
def generate_random_word(word_length, format='full', seed=0):
    '''
    Generate a random word with a specific length.
    ---------
    INPUTS:
    word_length: int; number of letters in the generated word.
    seed: int; random number generating seed; default 0.

    OUTPUTS:
    word: str; a word with letters in the alphabet and the specific length.
    '''
    if format not in ["full", "quad", "oct"]:
        print("Error, bad format!")
        raise ValueError
    word = ''
    for i in range(word_length):
        random.seed(seed + i)
        word += alphabet[int(len(alphabet) * random.random())]
    if format == 'full':
        return word
    if format == 'quad':
        print(word, quad_prefix)
        random.seed((seed + 100) * 10)  # must generate another random number
        prefix = quad_prefix[int(len(quad_prefix) * random.random())]
        return prefix + word
    if format == 'oct':  # not yet implemented
        # random.seed((seed+1000)*100) # must generate another random number
        # prefix = oct_prefix[int((len(quad_prefix) - 1) * random.random())]
        # return prefix+word
        return None

def get_rel_terms_in_symb(symb, fraction, rel, rel_slot='any', format='full', seed=0):
    '''
    Get the related term(s) in the given symbol according to the specified relation,
    for a fraction of words in the full symbol.
    ---------
    INPUTS:
    symb: dict.
    fraction: float; fraction of total words in the full symbol to be picked as samples;
              between 0 and 1, where 1 gets all words in the full symbol;
    rel: dict.
    rel_slot: str; one of the three choices 'first', 'final', 'any';
              if 'first' or 'final', the relation can only be at the first or final few slots of a word;
              if 'any', the relation can be at any slot of a word.
    seed: int; random number generating seed; default 0.

    OUTPUTS:
    rel_terms_list: list of dicts; each item in the list is a dict in the format of
                    {word: [symb_coeff, rel_coeff]}.
    '''
    if rel_slot not in {'first', 'initial', 'final', 'any'}: raise ValueError
    rel_terms_list_symb = []
    all_words = list(symb.keys())
    num_words_to_pick = int(len(all_words) * fraction)

    if num_words_to_pick <= 0:
        return []

    random.seed(seed)
    random_words = random.sample(all_words, num_words_to_pick)

    for word in random_words:
        rel_terms_list = get_rel_terms_in_symb_per_word(word, symb, rel, rel_slot=rel_slot, format=format)
        if rel_terms_list is None:
            return None
        for rel_terms in rel_terms_list:
            if (rel_terms) and (
                    not any([rel_terms == existing_rel_terms for existing_rel_terms in rel_terms_list_symb])):
                rel_terms_list_symb.append(rel_terms)

    return rel_terms_list_symb
def get_dihedral_terms_in_symb(word, symb, count_coeffs=False, failsymb=None):
    '''
    Get the {word: coeff} of all the dihedral images of a given word in a symbol;
    if word (or its images) is not in the symb, then coeff=0.
    ---------
    INPUTS:
    word: str.
    symb: dict.
    count_coeffs: bool; whether to output the number of occurances of a certain coeff among all dihedral images of a word;
                  default False.

    OUTPUTS:
    images_coeffs: dict; all terms dihedrally related to the given word in the symbol.
    unique_coeffs_counts: dict; number of occurances of each unique coeff among the six dihedral images of a word;
                          e.g., {'16' (coeff): 6 (# of occurances)}; return only if count_coeffs=True.
    '''
    images = get_dihedral_images(word)
    images_coeffs = {}
    for image in images:
        images_coeffs.update({image: get_coeff_from_word(image, symb)})
        if failsymb and image in failsymb: return None

    if not count_coeffs:
        return images_coeffs

    unique_coeffs, counts = np.unique(np.array(list(images_coeffs.values())), return_counts=True)
    return images_coeffs, dict(zip(unique_coeffs, counts))
def get_rel_instances_in_symb(rel_instance_list, symb):
    '''
    Given a list of relation instances, get their term(s) in the given symbol.
    First step for relation-level relation check.
    ---------
    INPUTS:
    rel_instance_list: list of dicts; outputs of functions `generate_rel_instances`.
    symb: dict; symbol at a given loop order; can be the ground truth or the model predictions.

    OUTPUTS:
    rel_instance_in_symb_list: list of dicts; each item in the list is a dict corresponding to one relation instance
                               in the input rel_instance_list; key is word and value is [symb_coeff, rel_coeff];
                               if the word is not in symb, then its symb_coeff is automatically 0.
    '''
    rel_instance_in_symb_list = []
    if rel_instance_list is None:
        return None

    for rel_instance in rel_instance_list:
        rel_instance_in_symb = {}
        for word, rel_coeff in rel_instance.items():
            symb_coeff = get_coeff_from_word(word, symb)
            rel_instance_in_symb.update({word: [symb_coeff, rel_coeff]})
        rel_instance_in_symb_list.append(rel_instance_in_symb)

    return rel_instance_in_symb_list
def update_rel_instances_in_symb(rel_instance_in_symb_list, symb):
    '''
    Given a list of relation instances in the format {word: [symb_coeff, rel_coeff]},
    update symb_coeff by the new input symbol.
    ---------
    INPUTS:
    rel_instance_list: list of dicts; format: {word: [symb_coeff, rel_coeff]}.
    symb: dict; symbol at a given loop order; can be the ground truth or the model predictions.

    OUTPUTS:
    rel_instance_in_symb_list: list of dicts; same as the output of function 'get_rel_instances_in_symb';
                               each item in the list is a dict corresponding to one relation instance
                               in the input rel_instance_list; key is word and value is [symb_coeff, rel_coeff].
    '''
    rel_instance_in_symb_list_updated = []
    if rel_instance_in_symb_list is None:
        return None

    for rel_instance in rel_instance_in_symb_list:
        rel_instance_in_symb = {}
        for word, coeff_pair in rel_instance.items():
            symb_coeff = get_coeff_from_word(word, symb)
            rel_instance_in_symb.update({word: [symb_coeff, coeff_pair[1]]})
        rel_instance_in_symb_list_updated.append(rel_instance_in_symb)

    return rel_instance_in_symb_list_updated
def replace_trivial0_terms(symb, return_symb=False):
    '''
    Find all the trivial zero words in a given symbol (assume valid and in full format)
    and manually force their coeffs to be zero, regardless of the original coeffs in the given symbol.
    ---------
    INPUTS:
    symb: dict.
    return_symb: bool; whether to return the full given symbol
                       with its trivial zero terms updated to have coeff=0;
                       default False.
    OUTPUTS:
    if return_symb == False, then
        trivial0_terms: dict; dict for trivial zero terms with {'word': 0} in the given symbol.
    if return_symb == True, then
        symb_updated: dict; full input symbol with its trivial zero terms updated to have coeff=0.
    '''
    trivial0_terms = {}
    symb_updated = symb.copy()
    for word in symb:
        if is_trivial0(word):
            trivial0_terms.update({word: 0})
            symb_updated.update({word: 0})

    if not return_symb:
        return trivial0_terms

    if return_symb:
        return symb_updated

def check_rel(rel_terms_list, return_rel_info=False, p_norm=None):
    '''
    Check if all relations in the input list of relations are satisfied, i.e., sum up to 0.
    Valid regardless of the sampling approach, i.e., word-oriented or relation-oriented.
    ---------
    INPUTS:
    rel_terms_list: list of dicts; format {word: [symb_coeff, rel_coeff]}.
    return_rel_info: bool; whether to return the detailed info of each relation,
                     including relation sum and number of non-trivial-zero terms in each relation;
                     default False.
    p_norm: float or None; p is the average accuracy of model prediction;
            goal is to normalize the accuracy by the number of terms in a relation; default None.

    OUTPUTS:
    percent: float; percent of correct relations, i.e., those that sum up to 0; between 0 and 1 (1: all correct).
    relsum_list: list; each item in the list is a int/float for one instance,
                 where the int/float corresponds to the rel sum of that particular instance (should all be 0 if correct);
                 return only if return_info_info=True.
    relnontrivial0_list: list; each item in the list is an int for one instance,
                 where the int is the number of non-trivial-zero words in that particular instance;
                 return only if return_rel_info=True.
    '''
    nterm = len(rel_terms_list[0])
    if rel_terms_list is None:
        return None

    relsum_list, relnontrivial0_list = zip(*[(elem) for elem in get_relsum_and_nzero(rel_terms_list, nterm, p_norm)])

    if not relsum_list:
        percent = None
    else:
        percent = relsum_list.count(0) / len(relsum_list)

        if p_norm:
            percent /= p_norm ** nterm

    if return_rel_info:
        return percent, relsum_list, relnontrivial0_list

    return percent
def check_coeffs_in_rel(rel_terms_list, symb_truth_list, return_counts=False, require_satisfied=True):
    '''
    Check if all relations in the input list of relations are satisfied, i.e., sum up to 0.
    Valid regardless of the sampling approach, i.e., word-oriented or relation-oriented.
    ---------
    INPUTS:
    rel_terms_list: list of dicts; format {word: [symb_coeff, rel_coeff]}.
    symb_truth: list of dicts;format {word: [symb_coeff, rel_coeff]}. the truth symbol against which the correctness of
                symb_coeff predicted by the model for each word is checked.
    return_counts: bool; whether to return the count of correct coeffs in each rel instance;
                    default False.

    OUTPUTS:
    percent_allcorrect: float; percent of correct relations with all its word coeffs also correct.
    percent_magcorrect: float; percent of correct relations with all its word coeffs magnitude correct.
    correct_coeffs_in_rel_list: list of lists; format [[bool, int, int], [], ...], where bool suggests if the current relation
                            instance is correct, the first int indicates how many of its word coeffs are right,
                            and the second int suggests how many word coeffs are magnitude correct;
                            return only if return_counts=True.
    '''
    correct_coeffs_in_rel_list = []
    n_allcorrect_rel, n_magcorrect_rel, n_signcorrect_rel = 0, 0, 0
    n_allcorrect_norel, n_magcorrect_norel, n_signcorrect_norel = 0, 0, 0
    if rel_terms_list is None:
        return None

    for rel_terms, symb_truth in zip(rel_terms_list, symb_truth_list):
        relsum = 0
        n_allcorrect, n_magcorrect, n_signcorrect = 0, 0, 0
        nterm = len(rel_terms)

        for (key, value), (truth_key, truth_value) in zip(rel_terms.items(), symb_truth.items()):
            if value[0] == None:  # invalid symb_coeff
                relsum = -1
            else:
                try:
                    relsum += np.prod(value)
                except ArithmeticError:
                    print("overflow error! setting rel sum to -1")
                    relsum = -1

            if value[0] != None:
                if value[0] == truth_value[0]:
                    n_allcorrect += 1

                if np.abs(value[0]) == np.abs(truth_value[0]):
                    n_magcorrect += 1

                if np.sign(value[0]) == np.sign(truth_value[0]):
                    n_signcorrect += 1

        if relsum == 0:
            rel_correct = True
        else:
            rel_correct = False

        if rel_correct == True and n_allcorrect == nterm:
            n_allcorrect_rel += 1

        if rel_correct == True and n_magcorrect == nterm:
            n_magcorrect_rel += 1

        if rel_correct == True and n_signcorrect == nterm:
            n_signcorrect_rel += 1

        if n_allcorrect == nterm:
            n_allcorrect_norel += 1

        if n_magcorrect == nterm:
            n_magcorrect_norel += 1

        if n_signcorrect == nterm:
            n_signcorrect_norel += 1

        correct_coeffs_in_rel_list.append([rel_correct, n_allcorrect, n_magcorrect, n_signcorrect])
    if not correct_coeffs_in_rel_list:
        percent_allcorrect, percent_magcorrect, percent_signcorrect = None, None, None
    else:
        if require_satisfied:
            percent_allcorrect = n_allcorrect_rel / len(correct_coeffs_in_rel_list)
            percent_magcorrect = n_magcorrect_rel / len(correct_coeffs_in_rel_list)
            percent_signcorrect = n_signcorrect_rel / len(correct_coeffs_in_rel_list)
        else:
            percent_allcorrect = n_allcorrect_norel / len(correct_coeffs_in_rel_list)
            percent_magcorrect = n_magcorrect_norel / len(correct_coeffs_in_rel_list)
            percent_signcorrect = n_signcorrect_norel / len(correct_coeffs_in_rel_list)
    if return_counts:
        return percent_allcorrect, percent_magcorrect, percent_signcorrect, correct_coeffs_in_rel_list

    return percent_allcorrect, percent_magcorrect, percent_signcorrect
def is_trivial0(word):
    '''
    Check if a given word (assuming valid and in full format) is a trivial zero word.
    ---------
    INPUTS:
    word: str.
    OUTPUTS:
    True/False: bool.
    '''
    for rel in first_entry_rel_table:  # prefix rule
        if word[0] in rel:
            return True

    for rel in final_entries_rel_table[:3]:  # suffix rule
        if word[-1] in rel:
            return True

    for rel in get_rel_table_dihedral(double_adjacency_rel_table):  # adjacency rule
        for rel_key in rel:
            if rel_key in word:
                return True

    return False
def get_rel_table_dihedral(rel_table):
    '''
    Given a relation table, output all the dihedral images of each relation,
    where duplicated relations are removed.
    ---------
    INPUTS:
    rel_table: list of dicts; e.g., one specific linear relations look-up table.

    OUTPUTS:
    unique_rel_table_dihedral: list of dicts; each item in the list is a dict corresponding to the dihedral images of one relation in the table;
                        total length of the list: nterms of the original table * 6 - duplicated terms.
    '''
    rel_table_dihedral = []
    unique_rel_table_dihedral = []
    for rel in rel_table:
        nterm = len(rel)  # number of terms in the relation
        term_list = list(rel.keys())

        term_list_dihedral = [get_dihedral_images(term) for term in term_list]
        for i in range(len(term_list_dihedral[0])):
            rel_dihedral = {}
            for iterm in range(nterm):
                rel_dihedral.update({term_list_dihedral[iterm][i]: rel[term_list[iterm]]})

            rel_table_dihedral.append(rel_dihedral)

    return dropdups(rel_table_dihedral)

    #seen_rel_table_dihedral = set()

    #    for d in rel_table_dihedral:
    #dict_tuple = tuple(sorted(d.items()))
    #        if dict_tuple not in seen_rel_table_dihedral:
    #            seen_rel_table_dihedral.add(dict_tuple)
    #            unique_rel_table_dihedral.append(d)

    #return unique_rel_table_dihedral
def get_rel_terms_in_symb_per_word(word, symb, rel, rel_slot='any', format='full'):
    '''
    Given a word, get the related term(s) in the given symbol according to the specified relation.
    This serves as an intermediate step to get related term(s) for more words in the full symbol.
    ---------
    INPUTS:
    word: str.
    symb: dict.
    rel: dict.
    rel_slot: str; one of the three choices 'first', 'final', 'any';
              if 'first' or 'final', the relation can only be at the first or final few slots of a word;
              if 'any', the relation can be at any slot of a word.

    OUTPUTS:
    rel_terms_list: list of dicts; each item in the list is a dict in the format of
                    {'bbbf' (word): [16 (coeff), 1 (rel coeff)]}.
    '''
    if rel_slot not in {'first', 'initial', 'final', 'any'}: raise ValueError
    rel_terms_list = []
    nterm = len(rel)
    nletter = len(next(iter(rel)))  # number of letters in each term
    if format not in ["full", "quad", "oct"]:
        print("Error, bad format!")
        raise ValueError

    # first entry relation
    if rel_slot == 'first':
        rel_terms = {}

        if format == 'full':
            if word[:nletter] not in rel:
                return rel_terms_list
            else:
                rel_terms.update({word: [get_coeff_from_word(word, symb), rel[word[:nletter]]]})
                rel_terms_list.append(rel_terms)
                return rel_terms_list
        else:  # compact formats
            if word[1:nletter] not in rel:  # ignore the prefix
                return rel_terms_list
            else:
                rel_terms.update({word: [get_coeff_from_word(word, symb), rel[word[1:nletter]]]})
                rel_terms_list.append(rel_terms)
                return rel_terms_list

    # first entry relation
    if rel_slot == 'initial':
        rel_terms = {}

        if format == 'full':
            if word[:nletter] not in rel:
                return rel_terms_list
            else:
                rel_terms.update({word: [get_coeff_from_word(word, symb), rel[word[:nletter]]]})
                rel_terms_list.append(rel_terms)
                return rel_terms_list
        else:  # compact formats
            if word[1:nletter] not in rel:  # ignore the prefix
                return rel_terms_list
            else:
                rel_terms.update({word: [get_coeff_from_word(word, symb), rel[word[1:nletter]]]})
                rel_terms_list.append(rel_terms)
                return rel_terms_list

    # final entries relation
    if rel_slot == 'final':
        rel_terms = {}

        if format == 'full':
            if word[-nletter:] not in rel:
                return rel_terms_list

            for key in rel:
                if word[-nletter:] == key:
                    pre_subword = word[:-nletter]
                    for key_rel in rel:
                        word_rel = pre_subword + key_rel
                        rel_terms.update({word_rel: [get_coeff_from_word(word_rel, symb), rel[key_rel]]})
                    rel_terms_list.append(rel_terms)
                    return rel_terms_list
        else:  # compact formats: no final entries relations
            return None

    # relations that can be at any position slot
    if rel_slot == 'any':

        key_rel_list = list(rel.keys())

        if format == 'full':
            word = word
            prefix = ''
        else:  # compact formats
            prefix = word[0]
            word = word[1:]

        if not any(key_rel in word for key_rel in key_rel_list):
            # rel_terms = {}
            # rel_terms_list.append(rel_terms)
            return rel_terms_list

        for key_rel in rel:
            if key_rel in word:
                start_pos_list = [i.start() for i in find_all(word, key_rel)]
                rel_terms_pos = []

                for start_pos in start_pos_list:
                    rel_terms = {}

                    pre_subword = word[:start_pos]
                    post_subword = word[start_pos + nletter:]

                    for key_rel in rel:
                        if format == 'full':
                            word_rel = pre_subword + key_rel + post_subword
                        else:  # compact formats
                            word_rel = prefix + pre_subword + key_rel + post_subword
                        rel_terms.update({word_rel: [get_coeff_from_word(word_rel, symb), rel[key_rel]]})

                    rel_terms_pos.append(rel_terms)

                rel_terms_list.append(rel_terms_pos)

        return list(itertools.chain(*rel_terms_list))

def get_relsum_and_nzero(rel_terms_list, nterm, p_norm):
    for rel_terms in rel_terms_list:
        relsum = 0
        n_nontrivial0_term = 0
        for key, value in rel_terms.items():
            if value[0] == None:  # invalid symb_coeff
                relsum = -1
            else:
                try:
                    relsum += np.prod(value)
                except ArithmeticError:
                    print("overflow error! setting rel sum to -1")
                    relsum = -1

            if not is_trivial0(key):
                n_nontrivial0_term += 1

        if p_norm:
            try:
                relsum /= p_norm ** nterm
            except ArithmeticError:
                print("overflow error! setting rel sum to -1")
                relsum = -1
        yield relsum, n_nontrivial0_term
def gen_let(inletter,type="next"):
    # hardcode adjacency rules to generate nontrivial zeroes
    assert type in {"next","last","first"}
    tmat={'a':{"first":['a', 'b', 'c'],
        "next":['a', 'b', 'c', 'e', 'f'],
        "last":['e','f']},
    'b': {"first":['a', 'b', 'c'],
        "next":['a', 'b', 'c', 'd', 'f'],
        "last":['d','f']},
    'c':{"first":['a','b','c'],
        "next":['a', 'b', 'c', 'd', 'e'],
        "last":['d','e']},
    'd':{"first":['b','c'],
        "next":['b', 'c', 'd'],
        "last":['e']},
    'e':{"first":['a','c'],
        "next":['a', 'c', 'e'],
        "last":['e']},
    'f':{"first":['a','b'],
        "next":['a', 'b', 'f'],
        "last":['f']}
    }
    letter = tmat[inletter][type][int(len(tmat[inletter][type]) * random.random())]
    return letter
def gen_valid_substr(to_gen, input=None, suffix=False):
    # generate a valid substring. If input is given, build a string that is compatible with it.
    # If 'suffix', gen a valid substring that can be reversed and prepended to input
    # otherwise, gen a valid substring that can be appended to input
    letter = ''
    if input:
        if suffix:
            letter = input[0]
        else:
            letter = input[-1]
    for i in range(to_gen):
        if not suffix:
            if i == 0:
                mylist = ['a', 'b', 'c']
                letter = mylist[int(len(mylist) * random.random())]
        if suffix and (i == to_gen - 1):
            letter = gen_let(letter,"first")
        else:
            letter = gen_let(letter,"next")
        yield letter
