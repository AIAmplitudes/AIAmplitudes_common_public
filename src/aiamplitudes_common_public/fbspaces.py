"""
Front-space and back-space basis operations.

Front/back spaces decompose symbol prefixes/suffixes into independent basis
elements using the physics relations. Two representations are supported:

  - Permissive (P): Overcomplete basis of SB(...) linear combinations.
    Elements named FP_w_i / BP_w_i. Loaded from frontspace/backspace files.
  - Restrictive (R): Minimal independent basis of E(...) or SB(...) elements.
    Elements named FR_w_i / BR_w_i. Loaded from multiinitial/multifinal files.

Each basis has a forward dict (name -> letter expansion) and a flip dict
(letter string -> name) for bidirectional lookup.

Dimensions are read from the selected basis and cached on first load.
F_number[w] and B_number[w] expose the restrictive dimensions;
F_perm_number[w] and B_perm_number[w] expose the permissive dimensions.
NF_rels[w] and NB_rels[w] similarly expose cached restrictive relation counts.

Also provides symbol expansion/compression and relation lookup functions
that are shared by both common_public and common_dev.
"""

import re
import os
from fractions import Fraction
from numbers import Integral
from collections import Counter
from aiamplitudes_common_public.commonclasses import Symb
from aiamplitudes_common_public.file_readers import readSymb, readFile, SB_to_dict, SB_to_dict_fraction
from aiamplitudes_common_public.download_data import relpath

# Shared process-local caches. Keys include the source file, so different bases
# at the same weight never share a dimension by accident.
_basis_cache = {}
basis_sizes = {}
_relation_cache = {}
relation_counts = {}
_basis_cache_clearers = []


class _BasisDimensions:
    """Weight-indexed view of the selected basis's cached dimensions."""

    def __init__(self, space, rp):
        self.space = space
        self.rp = rp
        self._sizes = {}

    def __getitem__(self, weight):
        if not isinstance(weight, Integral) or isinstance(weight, bool) or weight < 0:
            raise ValueError("weight must be a nonnegative integer")
        if weight not in self._sizes:
            from aiamplitudes_common_public import basis_size
            self._sizes[weight] = basis_size(weight, self.space, self.rp)
        return self._sizes[weight]


F_number = _BasisDimensions('front', 'R')
B_number = _BasisDimensions('back', 'R')
F_perm_number = _BasisDimensions('front', 'P')
B_perm_number = _BasisDimensions('back', 'P')


class _RelationCounts:
    """Weight-indexed view of the selected restrictive relation counts."""

    def __init__(self, space):
        self.space = space
        self._counts = {}

    def __getitem__(self, weight):
        if not isinstance(weight, Integral) or isinstance(weight, bool) or weight < 0:
            raise ValueError("weight must be a nonnegative integer")
        if weight not in self._counts:
            from aiamplitudes_common_public import relation_count
            self._counts[weight] = relation_count(weight, self.space)
        return self._counts[weight]


NF_rels = _RelationCounts('front')
NB_rels = _RelationCounts('back')


def clear_basis_cache():
    """Forget bases, relations, their counts, and derived data after file changes."""
    _basis_cache.clear()
    basis_sizes.clear()
    _relation_cache.clear()
    relation_counts.clear()
    for counts in (NF_rels, NB_rels):
        counts._counts.clear()
    for dimensions in (F_number, B_number, F_perm_number, B_perm_number):
        dimensions._sizes.clear()
    _header_cache.clear()
    _preload_cache.clear()
    for clear in _basis_cache_clearers:
        clear()


def _basis_filepath(prefix):
    return os.path.abspath(os.path.join(relpath, os.fspath(prefix)))

_header_cache = {}

def _parse_file_header(filepath, kind):
    """Parse the 'init' or 'fin' relations and basis names from a file header.

    Returns (relnames, spacenames) dicts keyed by weight (1-indexed).
    """
    cache_key = (filepath, kind)
    if cache_key in _header_cache:
        return _header_cache[cache_key]

    relnames = {}
    spacenames = {}

    with open(filepath, 'rt') as f:
        full = f.read()

    # Parse fin_rels / init_rels: list of label names
    m = re.search(rf'{kind}_rels\s*:=\s*\[([^\]]+)\]', full)
    if m:
        rels_list = [s.strip() for s in m.group(1).split(',') if s.strip()]
        for i, name in enumerate(rels_list):
            relnames[i + 1] = name

    # Parse fin_list / init_list entries: fin_list[w] := labelname :
    # Skip commented-out lines (preceded by # on the same line)
    for m in re.finditer(rf'{kind}_list\[(\d+)\]\s*:=\s*(\w+)', full):
        line_start = full.rfind('\n', 0, m.start()) + 1
        if full[line_start:m.start()].lstrip().startswith('#'):
            continue
        w = int(m.group(1))
        spacenames[w] = m.group(2)

    # Older multiinitial/multifinal files omit some or all metadata headers.
    # Discover their canonical section labels; numeric suffixes are names, not
    # trusted dimensions/counts. Explicit metadata above remains authoritative.
    stems = ("single", "double", "triple", "quad", "quint", "hex", "hept",
             "oct", "non", "dec")
    prefix = "i" if kind == "init" else ""
    section = re.compile(
        rf"^\s*{prefix}({'|'.join(stems)})(indep|rels)\d+\s*:=", re.MULTILINE
    )
    for match in section.finditer(full):
        weight = stems.index(match.group(1)) + 1
        label = match.group(0).split(":=", 1)[0].strip()
        names = spacenames if match.group(2) == "indep" else relnames
        if weight not in names:
            names[weight] = label

    # Only keep relnames for weights that have a corresponding spacenames entry
    relnames = {w: v for w, v in relnames.items() if w in spacenames}

    _header_cache[cache_key] = (relnames, spacenames)
    return relnames, spacenames

def _get_brelnames(prefix):
    relnames, _ = _parse_file_header(_basis_filepath(prefix), 'fin')
    return relnames

def _get_bspacenames(prefix):
    _, spacenames = _parse_file_header(_basis_filepath(prefix), 'fin')
    return spacenames

def _get_frelnames(prefix):
    relnames, _ = _parse_file_header(_basis_filepath(prefix), 'init')
    return relnames

def _get_fspacenames(prefix):
    _, spacenames = _parse_file_header(_basis_filepath(prefix), 'init')
    return spacenames

# Legacy hardcoded dicts (kept for backward compat with w=7 front-space split)
frelnames = {1: 'isinglerels3',
             2: 'idoublerels9',
             3: 'itriplerels33',
             4: 'iquadrels78',
             5: 'iquintrels180',
             6: 'ihexrels406',
             7: {'zero':'iheptrels_z', 'nonzero':'iheptrels_c'},
             8: 'ioctrels2079',
            }


def get_perm_fspace(w, prefix='frontspace'):
    """Load the permissive front-space basis at weight w.

    Returns (basedict, flipdict) where:
      basedict: {'FP_w_i': {letter_string: coefficient, ...}, ...}
      flipdict: {letter_string: {'FP_w_i': coefficient, ...}, ...} (reverse lookup)
    """
    assert os.path.isfile(_basis_filepath(prefix))
    mystr = ''.join(str.split(readSymb(_basis_filepath(prefix),'frontspace',w)))
    newstr = re.split(r":=|\[|\]", mystr)[4]
    dev = [elem + ")" if elem[-1] != ")" else elem for elem in newstr.split("),") if elem]
    basedict = {f'FP_{w}_{i+1}': SB_to_dict_fraction(el) for i, el in enumerate(dev)}
    flipdict = {}
    for elem, elemdict in basedict.items():
        for term, coef in elemdict.items():
            if term not in flipdict: flipdict[term] = {}
            flipdict[term][elem] = basedict[elem][term]
    return basedict, flipdict

def get_perm_bspace(w, prefix = 'backspace'):
    """Load the permissive back-space basis at weight w. Returns (basedict, flipdict)."""
    assert os.path.isfile(_basis_filepath(prefix))
    mystr = ''.join(str.split(readSymb(_basis_filepath(prefix), 'backspace', w)))
    newstr = re.split(r":=|\[|\]", mystr)[4]
    dev = [elem + ")" if elem[-1] != ")" else elem for elem in newstr.split("),") if elem]
    basedict = {f'BP_{w}_{i+1}': SB_to_dict_fraction(el) for i, el in enumerate(dev)}
    flipdict = {}
    for elem, elemdict in basedict.items():
        for term, coef in elemdict.items():
            if term not in flipdict: flipdict[term] = {}
            flipdict[term][elem] = basedict[elem][term]
    return basedict, flipdict

def get_rest_bspace(w, prefix = 'phi2multifinal_E'):
    """Load the restrictive back-space basis at weight w from the multifinal file.

    Returns (flip, myd) where:
      flip: {'BR_w_i': letter_string, ...} (basis name -> letters)
      myd: {letter_string: 'BR_w_i', ...} (letters -> basis name)
    """
    assert os.path.isfile(_basis_filepath(prefix))
    spacenames = _get_bspacenames(prefix)
    if w not in spacenames:
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(spacenames.keys())})")
    res=readSymb(_basis_filepath(prefix),str(spacenames[w]))
    if not res.strip():
        raise ValueError(
            f"Weight {w} is declared in '{prefix}' (label '{spacenames[w]}') but no matching "
            f"data block was found in the file body -- it may be stored under a different or "
            f"compressed label (e.g. a '_G'-suffixed block)."
        )
    myindeps = [elem for elem in re.split(r":=\[|E\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem]
    myd = {elem: f'BR_{w}_{i+1}' for i, elem in enumerate(myindeps)}
    flip = {f'BR_{w}_{i+1}': elem for i, elem in enumerate(myindeps)}
    return flip, myd

def get_rest_fspace(w, prefix='multiinitial_E'):
    """Load the restrictive front-space basis at weight w from the multiinitial file.

    Returns (flip, myd) with the same structure as get_rest_bspace but using FR_ names.
    """
    assert os.path.isfile(_basis_filepath(prefix))
    spacenames = _get_fspacenames(prefix)
    if w not in spacenames:
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(spacenames.keys())})")
    res=readSymb(_basis_filepath(prefix),str(spacenames[w]))
    if not res.strip():
        raise ValueError(
            f"Weight {w} is declared in '{prefix}' (label '{spacenames[w]}') but no matching "
            f"data block was found in the file body -- it may be stored under a different or "
            f"compressed label (e.g. a '_G'-suffixed block)."
        )
    myindeps = [elem for elem in re.split(r":=\[|SB\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem]
    myd = {elem: f'FR_{w}_{i+1}' for i, elem in enumerate(myindeps)}
    flip = {f'FR_{w}_{i+1}': elem for i, elem in enumerate(myindeps)}
    return flip, myd

def getBrel_eqs(f, w, prefix='phi2multifinal_E'):
    """Extract back-space relation equations at weight w from an open file handle."""
    relnames = _get_brelnames(prefix)
    if w not in relnames:
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(relnames.keys())})")
    res = readFile(f, str(relnames[w]) + ' :=')
    if not res.strip():
        raise ValueError(
            f"Weight {w} is declared in '{prefix}' (label '{relnames[w]}') but no matching "
            f"relation block was found in the file body -- it may be stored under a different "
            f"or compressed label."
        )
    out = _parse_rel_block(res)
    out = _resolve_ops(out, f)
    return out

def _parse_rel_block(res):
    """Parse a single relation block string into a list of equation strings."""
    # Remove commas only inside parentheses (SB(a,b,c) -> SB(abc), E(a,b) -> E(ab))
    def _strip_inner_commas(s):
        out = []
        depth = 0
        for ch in s:
            if ch == '(':
                depth += 1
                out.append(ch)
            elif ch == ')':
                depth -= 1
                out.append(ch)
            elif ch == ',' and depth > 0:
                pass
            else:
                out.append(ch)
        return ''.join(out)

    cleaned = _strip_inner_commas(res)
    return [re.sub(r'\s+', '', elem) for elem in
            re.split(r":=\s*\[|,|\]\s*:", cleaned)[1:] if elem.strip()]

def _resolve_ops(elems, f):
    """Resolve any op(...) references by reading the referenced sections from f."""
    out = []
    for elem in elems:
        m = re.match(r'^op\((\w+)\)$', elem)
        if m:
            ref_name = m.group(1)
            f.seek(0)
            ref_res = readFile(f, ref_name + ' :=')
            if not ref_res.strip():
                raise ValueError(f"op({ref_name}) reference could not be resolved -- "
                                  f"no block labeled '{ref_name}' found in the file body.")
            out.extend(_parse_rel_block(ref_res))
        else:
            out.append(elem)
    return out

def getFrel_eqs(f, w, prefix='multiinitial_E'):
    """Extract front-space relation equations at weight w. Handles split files via op()."""
    relnames = _get_frelnames(prefix)
    if w not in relnames and not isinstance(frelnames.get(w), dict):
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(relnames.keys())})")
    # Honor the selected file's declaration, including custom weight-7 bases.
    # Older files may only contain the split zero/nonzero blocks below.
    if w in relnames:
        f.seek(0)
        res = readFile(f, str(relnames[w]) + ' :=')
        if res.strip():
            return _resolve_ops(_parse_rel_block(res), f)
    if isinstance(frelnames.get(w), dict):
        zeros = str(frelnames[w]['zero'])
        nonzeros = str(frelnames[w]['nonzero'])
        f.seek(0)
        res_zero = readFile(f, zeros + ' :=')
        f.seek(0)
        res_nonzero = readFile(f, nonzeros + ' :=')
        if not res_zero.strip() or not res_nonzero.strip():
            raise ValueError(
                f"Weight {w} split relation blocks ('{zeros}'/'{nonzeros}') not both found in "
                f"'{prefix}' -- one or both blocks are missing from the file body."
            )
        out_zero = _parse_rel_block(res_zero)
        out_nonzero = _parse_rel_block(res_nonzero)
        out = _resolve_ops(out_zero + out_nonzero, f)
    else:
        raise ValueError(
            f"Weight {w} is declared in '{prefix}' (label '{relnames[w]}') but no matching "
            f"relation block was found in the file body -- it may be stored under a "
            f"different or compressed label."
        )
    return out

def rel_to_dict(relstring, bspace=True):
    """Parse a relation string into a nested dict.

    e.g. 'E(abc)=-2*E(def)+4*E(bcd)' -> {'abc': {'def': -2, 'bcd': 4}}
    The outer key is the dependent element; inner dict gives the independent
    elements and their coefficients in the expansion.
    """
    if not relstring or relstring.strip() in ('', 'NULL'):
        return {None: None}

    def expandcoef(c): return Fraction(c + '1') if (len(c) == 1 and not c.isnumeric()) else Fraction(c)

    if bspace:
        newstring = [elem for elem in re.split(r"=|E\(|\)", re.sub('#', '',
                                                                  re.sub('[,*]', '', relstring))) if elem]
    else:
        newstring = [elem for elem in re.split(r"=|SB\(|\)", re.sub('#', '',
                                                                  re.sub('[,*]', '', relstring))) if elem]

    if len(newstring) == 0: return {None: None}
    if len(newstring) == 1: return {newstring[0]: {None: 0}}
    if newstring[1] == '0':
        reldict = {None: 0}
    else:
        if newstring[1].isalpha():
            eq = ['+'] + newstring[1:]
        else:
            eq = newstring[1:]

        if len(eq) == 1:
            reldict = {eq[0]: 1}
        else:
            reldict = {k: expandcoef(c) for i, (c, k) in enumerate(zip(eq, eq[1:])) if i % 2 == 0}

    return {newstring[0]: reldict}

def _relation_key(w, space, directory, prefix):
    if not isinstance(w, Integral) or isinstance(w, bool) or w < 0:
        raise ValueError("weight must be a nonnegative integer")
    if space not in {'front', 'back'}:
        raise ValueError("space must be 'front' or 'back'")
    filepath = os.path.abspath(os.path.join(directory, os.fspath(prefix)))
    return space, filepath, int(w)


def _load_relations(key):
    if key not in _relation_cache:
        space, filepath, w = key
        if w == 0:
            relations = {}
        else:
            reader = getFrel_eqs if space == 'front' else getBrel_eqs
            with open(filepath, 'rt') as stream:
                relations = {
                    word: expansion
                    for equation in reader(stream, w, prefix=filepath)
                    for word, expansion in rel_to_dict(equation, space == 'back').items()
                    if word
                }
        relation_counts[key] = len(relations)
        _relation_cache[key] = relations
    return _relation_cache[key]


def _relation_count(w, space, directory, prefix):
    key = _relation_key(w, space, directory, prefix)
    if key not in relation_counts:
        _load_relations(key)
    return relation_counts[key]


def get_brels(w, relpath, prefix='phi2multifinal_E'):
    """Load cached back relations and record their count. Treat the dict as read-only."""
    return _load_relations(_relation_key(w, 'back', relpath, prefix))


def get_frels(w, relpath, prefix='multiinitial_E'):
    """Load cached front relations and record their count. Treat the dict as read-only."""
    return _load_relations(_relation_key(w, 'front', relpath, prefix))


def all_perm_bspaces(relpath):
    perm_bspace, perm_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        perm_bspace |= get_perm_bspace(relpath, i)

def all_rest_bspaces(relpath):
    rest_bspace, rest_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        rest_bspace |= get_rest_bspace(relpath, i)


########################################################
# Functions to expand and compress symbols between F/B-space and letter forms
########################################################

def expand_elem(key, opt="all", loaded=None):
    """Expand a compressed key (with F/B-space names) to full letter-basis Symb.

    Each basis element maps to a linear combination of letter strings, so the
    result is a Symb of {full_letter_key: product_of_coefficients}.
    """
    def proc_elem(myelem):
        for entry in myelem.split('@'):
            yield entry
    elem = [out for out in proc_elem(key)]
    if opt == "front":
        #fweight, bweight = int(elem[0].split('_')[1]), None
        myfspace=loaded["fspace"]
        fkey = elem[0]
        if fkey not in myfspace:
            parts = fkey.split('_')
            fkey = f'FPD_{parts[1]}_{parts[2]}'
        mid = ''.join(elem[1:])
        exp_dict = {f'{k}{mid}':v for k, v in myfspace[fkey].items()}
    elif opt == "back":
        #fweight, bweight = None, int(elem[-1].split('_')[1])
        mybspace = loaded["bspace"]
        mid=''.join(elem[:-1])
        bkey = elem[-1]
        if bkey not in mybspace:
            parts = bkey.split('_')
            bkey = f'BPD_{parts[1]}_{parts[2]}'
        exp_dict = {f'{mid}{k}':v for k, v in mybspace[bkey].items()}
    elif opt == "all":
        #fweight, bweight = int(elem[0].split('_')[1]), int(elem[-1].split('_')[1])
        myfspace = loaded["fspace"]
        mybspace = loaded["bspace"]
        mid=''.join(elem[1:-1])
        #print(elem, myfspace[elem[0]], mybspace[elem[-1]])
        exp_dict = {f'{k1}{mid}{k2}': v1*v2
                    for k1, v1 in myfspace[elem[0]].items()
                    for k2, v2 in mybspace[elem[-1]].items()}
    else:
        raise ValueError
    return Symb(exp_dict)

def _load_fspace(fForm, fweight):
    """Load the selected front basis and its associated relations."""
    from aiamplitudes_common_public import fspace, fspace_flip, fr_rels, basis_size
    loaded = {
        "fspace": fspace(fweight, fForm),
        "fspace_flip": fspace_flip(fweight, fForm),
        "fspace_size": basis_size(fweight, 'front', fForm),
    }
    if fForm == 'R' and fweight > 0:
        loaded["f_rels"] = fr_rels(fweight)
    return loaded


def _load_bspace(bForm, bweight, optrace='phi2'):
    """Load the selected back basis and its associated relations."""
    from aiamplitudes_common_public import bspace, bspace_flip, br_rels, basis_size
    loaded = {
        "bspace": bspace(bweight, bForm, optrace=optrace),
        "bspace_flip": bspace_flip(bweight, bForm, optrace=optrace),
        "bspace_size": basis_size(bweight, 'back', bForm, optrace=optrace),
    }
    if bForm == 'R' and bweight > 0:
        loaded["b_rels"] = br_rels(bweight, optrace=optrace)
    return loaded


_preload_cache = {}

def preload_fbspaces(seam, fForm, fweight, bForm, bweight, optrace='phi2'):
    """Pre-load F/B-space basis dicts and relations into a single dict for batch processing.

    Results are cached for reuse across calls with the same arguments.
    """
    cache_key = (seam, fForm, fweight, bForm, bweight, optrace)
    if cache_key in _preload_cache:
        return _preload_cache[cache_key]

    loaded = dict()

    if seam == 'front':
        loaded.update(_load_fspace(fForm, fweight))
    elif seam == 'back':
        loaded.update(_load_bspace(bForm, bweight, optrace=optrace))
    elif seam == 'all':
        loaded.update(_load_fspace(fForm, fweight))
        loaded.update(_load_bspace(bForm, bweight, optrace=optrace))
    else:
        raise ValueError

    _preload_cache[cache_key] = loaded
    return loaded


def get_elem_compression(pattern):
    """Detect the compression scheme from a single key string.

    Returns (fForm, fweight, seamsize, bForm, bweight).
    """
    parts = pattern.split('@')

    if len(parts) == 1:
        return None, None, len(parts[0]), None, None

    if '_' in parts[0] and '_' in parts[-1]:
        if len(parts) == 2: seamsize = 0
        else: seamsize= len(parts[1])
    elif '_' in parts[0]: seamsize=len(parts[1])
    elif '_' in parts[1]: seamsize = len(parts[0])
    else: seamsize = len(parts[0])

    first = parts[0].split('_')
    last = parts[-1].split('_')
    fweight, bweight = None, None
    fForm, bForm = None, None
    if 'P' in first[0] or 'R' in first[0]: fweight = int(first[1])
    if 'P' in last[0] or 'R' in last[0]: bweight = int(last[1])

    if 'PD' in first[0]:
        fForm = 'PD'
    elif 'P' in first[0]:
        fForm = 'P'
    elif 'R' in first[0]:
        fForm = 'PD'

    if 'PD' in last[0]:
        bForm = 'PD'
    elif 'P' in last[0]:
        bForm = 'P'
    elif 'R' in last[0]:
        bForm = 'PD'

    return fForm, fweight, seamsize, bForm, bweight

def get_compression(mySymb):
    """Detect the compression scheme of a Symb by examining its first key."""
    pattern=next(iter(mySymb))
    return get_elem_compression(pattern)

def expand_rform_elem(key, opt="all", loaded=None):
    """Expand a restrictive-form key into individual letter components (flat list).

    Returns (expanded_list, (fweight, bweight)).
    """
    def proc_elem(myelem):
        for entry in myelem.split('@'):
            if entry.isalpha():
                for i in entry:
                    yield i
            else:
                yield entry

    elem = [out for out in proc_elem(key)]
    if opt == "front":
        fweight, bweight = int(elem[0].split('_')[1]), None
        myfspace=loaded["fspace"] if loaded else get_rest_fspace(fweight)[0]
        expanded = [i for i in myfspace[elem[0]]] + elem[1:]
    elif opt == "back":
        fweight, bweight = None, int(elem[-1].split('_')[1])
        mybspace = loaded["bspace"] if loaded else get_rest_bspace(bweight)[0]
        expanded = elem[:-1] + [i for i in mybspace[elem[-1]]]
    elif opt == "all":
        fweight, bweight = int(elem[0].split('_')[1]), int(elem[-1].split('_')[1])
        myfspace = loaded["fspace"] if loaded else get_rest_fspace(fweight)[0]
        mybspace = loaded["bspace"] if loaded else get_rest_bspace(bweight)[0]

        expanded = [i for i in myfspace[elem[0]]] + elem[1:-1] + [i for i in mybspace[elem[-1]]]
    else:
        raise ValueError
    return expanded, (fweight, bweight)


def compress_rform_elem(key, fweight,bweight, seam,loaded=None):
    """Compress a letter list back to restrictive-form notation using flip dicts."""
    if seam in {"front", "all"}:
        fseam = fweight
        myfflip = loaded["fspace_flip"] if loaded else get_rest_fspace(fweight)[1]
        myf = [myfflip[''.join(key[:fseam])]]
    if seam in {"back", "all"}:
        bseam = len(key) - bweight
        mybflip = loaded["bspace_flip"] if loaded else get_rest_bspace(bweight)[1]
        myb = [mybflip[''.join(key[bseam:])]]

    if seam == "front":
        suffix_parts = key[fseam:]
        if suffix_parts and '_' in str(suffix_parts[-1]):
            l = myf + ([''.join(suffix_parts[:-1])] if len(suffix_parts) > 1 else []) + [suffix_parts[-1]]
        else:
            l = myf + [''.join(suffix_parts)] if suffix_parts else myf
    elif seam == "back":
        prefix_parts = key[:bseam]
        if prefix_parts and '_' in str(prefix_parts[0]):
            l = [prefix_parts[0]] + ([''.join(prefix_parts[1:])] if len(prefix_parts) > 1 else []) + myb
        else:
            l = [''.join(prefix_parts)] + myb
    elif seam == "all":
        l = myf + [''.join(key[fseam:bseam])] + myb
    return '@'.join([p for p in l if p])

def get_related(rel, seed, slot):
    """Given a seed word and relation, get all related words at the given slot."""
    seedword = list(seed)
    return {','.join(seedword[:slot] + [i for i in elem] + seedword[slot + len(elem):]): relcoef
            for elem, relcoef in rel.items()}


def get_as_indepsum(fullword, fweight=None, bweight=None, seam=None,
                    loaded=None):
    """Express a non-independent term as a sum of basis (independent) terms."""
    if isinstance(fullword, list):
        fullword = ''.join(fullword)
    if seam == "back":
        mybflip = loaded["bspace_flip"] if loaded else get_rest_bspace(bweight)[1]
        mybrels = loaded["b_rels"] if loaded else get_brels(bweight,relpath)
        body, belem = fullword[:-bweight], ''.join(fullword[-bweight:])
        if belem in mybflip:
            reldict = {compress_rform_elem(fullword, fweight, bweight, seam, loaded): 1}
        elif belem in mybrels:
            reldict = {compress_rform_elem(body + ''.join([i for i in k]), fweight, bweight, seam, loaded): v
                       for k, v in mybrels[belem].items() if k}
        else: reldict = {}

    elif seam == "front":
        myfflip = loaded["fspace_flip"] if loaded else get_rest_fspace(fweight)[1]
        myfrels = loaded["f_rels"] if loaded else get_frels(fweight,relpath)
        felem, body = ''.join(fullword[:fweight]), fullword[fweight:]
        if felem in myfflip:
            reldict = {compress_rform_elem(fullword, fweight, bweight, seam,loaded): 1}
        elif felem in myfrels:
            reldict = {compress_rform_elem(''.join([i for i in k]) + body, fweight, bweight, seam, loaded): v
                       for k, v in myfrels[felem].items() if k}
        else: reldict = {}
    else:
        print("Bad seam!")
        raise ValueError
    return (reldict)



def get_perm_fspace_wt6(filename='wt6_242_indep_symbols'):
    """Read the 242 independent weight-6 front-space symbols.

    The file contains a Maple-style list of SB(...) linear combinations with
    potentially fractional coefficients.

    Returns:
        List of 242 dicts, each mapping letter-strings to Fraction coefficients.
    """
    filepath = _basis_filepath(filename)
    with open(filepath, 'rt') as f:
        content = f.read()

    # Strip comments and extract the list between := [ ... ]:
    content = re.sub(r'#[^\n]*', '', content)
    match = re.search(r':=\s*\[(.*)\]\s*:', content, re.DOTALL)
    if not match:
        raise ValueError(f"Could not parse list from {filepath}")
    body = match.group(1)

    # Remove Maple line continuations (backslash + newline) then remaining newlines
    body = body.replace('\\\n', '').replace('\n', '')

    # Split on commas that separate top-level elements.
    # Each element is a sum of coeff*SB(...) terms.
    # Split carefully: commas inside SB() are not separators.
    elements = []
    depth = 0
    current = []
    for ch in body:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == ',' and depth == 0:
            elements.append(''.join(current).strip())
            current = []
        else:
            current.append(ch)
    if current:
        elements.append(''.join(current).strip())

    # Parse each element using SB_to_dict_fraction
    result = {}
    idx = 0
    for elem in elements:
        if not elem:
            continue
        # Add closing paren if stripped during split
        if elem.count('(') > elem.count(')'):
            elem += ')'
        d = SB_to_dict_fraction(elem)
        if d:
            idx += 1
            result[f'FP_6_{idx}'] = d

    return result
