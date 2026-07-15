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

Dimension tables (indexed by weight 0-8):
  B_number (back):  [1, 3, 6, 12, 24, 45, 85, 155, 279]
  F_number (front): [1, 3, 9, 21, 48, 108, 246, 555, 1251]

Also provides symbol expansion/compression and relation lookup functions
that are shared by both common_public and common_dev.
"""

import re
import os
from fractions import Fraction
from collections import Counter
from aiamplitudes_common_public.commonclasses import Symb
from aiamplitudes_common_public.file_readers import readSymb, readFile, SB_to_dict, SB_to_dict_fraction
from aiamplitudes_common_public.download_data import relpath

B_number= [1, 3, 6, 12, 24, 45, 85, 155, 279, None ] # dim of back space at each weight
F_number= [1, 3, 9, 21, 48, 108, 246, 555, 1251, None ] # dim of front space at each weight

NB_rels= [0, 3, 12, 24, 48, 99, 185, 355, 651, None ] # number of back-space relations
NF_rels= [0, 3, 9, 33, 78, 180, 402, 921, 2079, None] # number of front-space relations

_header_cache = {}

def _parse_file_header(filepath):
    """Parse fin_rels/init_rels and fin_list/init_list from a multifinal/multiinitial file header.

    Returns (relnames, spacenames) dicts keyed by weight (1-indexed).
    """
    if filepath in _header_cache:
        return _header_cache[filepath]

    relnames = {}
    spacenames = {}

    with open(filepath, 'rt') as f:
        full = f.read()

    # Parse fin_rels / init_rels: list of label names
    m = re.search(r'(?:fin_rels|init_rels)\s*:=\s*\[([^\]]+)\]', full)
    if m:
        rels_list = [s.strip() for s in m.group(1).split(',') if s.strip()]
        for i, name in enumerate(rels_list):
            relnames[i + 1] = name

    # Parse fin_list / init_list entries: fin_list[w] := labelname :
    # Skip commented-out lines (preceded by # on the same line)
    for m in re.finditer(r'(?:fin_list|init_list)\[(\d+)\]\s*:=\s*(\w+)', full):
        line_start = full.rfind('\n', 0, m.start()) + 1
        if full[line_start:m.start()].lstrip().startswith('#'):
            continue
        w = int(m.group(1))
        spacenames[w] = m.group(2)

    # Only keep relnames for weights that have a corresponding spacenames entry
    relnames = {w: v for w, v in relnames.items() if w in spacenames}

    _header_cache[filepath] = (relnames, spacenames)
    return relnames, spacenames

def _get_brelnames(prefix):
    relnames, _ = _parse_file_header(f'{relpath}/{prefix}')
    return relnames

def _get_bspacenames(prefix):
    _, spacenames = _parse_file_header(f'{relpath}/{prefix}')
    return spacenames

def _get_frelnames(prefix):
    relnames, _ = _parse_file_header(f'{relpath}/{prefix}')
    return relnames

def _get_fspacenames(prefix):
    _, spacenames = _parse_file_header(f'{relpath}/{prefix}')
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
    assert os.path.isfile(f'{relpath}/{prefix}')
    mystr = ''.join(str.split(readSymb(f'{relpath}/{prefix}','frontspace',w)))
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
    assert os.path.isfile(f'{relpath}/{prefix}')
    mystr = ''.join(str.split(readSymb(f'{relpath}/{prefix}', 'backspace', w)))
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
    assert os.path.isfile(f'{relpath}/{prefix}')
    spacenames = _get_bspacenames(prefix)
    if w not in spacenames:
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(spacenames.keys())})")
    res=readSymb(f'{relpath}/{prefix}',str(spacenames[w]))
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
    assert os.path.isfile(f'{relpath}/{prefix}')
    spacenames = _get_fspacenames(prefix)
    if w not in spacenames:
        raise KeyError(f"Weight {w} not available in '{prefix}' (max: {max(spacenames.keys())})")
    res=readSymb(f'{relpath}/{prefix}',str(spacenames[w]))
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
    # Override for known split case (label in header doesn't exist as a data block)
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
        res = readFile(f, str(relnames[w]) + ' :=')
        if not res.strip():
            raise ValueError(
                f"Weight {w} is declared in '{prefix}' (label '{relnames[w]}') but no matching "
                f"relation block was found in the file body -- it may be stored under a "
                f"different or compressed label."
            )
        out = _parse_rel_block(res)
        out = _resolve_ops(out, f)
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

def get_brels(w,relpath, prefix='phi2multifinal_E'):
    """Load all back-space relations at weight w as a flat {dependent: {indep: coef}} dict."""
    assert w > 0
    with open(f'{relpath}/{prefix}', 'rt') as f:
        return {k: v for j in getBrel_eqs(f, w, prefix=prefix) for k, v in rel_to_dict(j).items() if k}

def get_frels(w,relpath, prefix='multiinitial_E'):
    """Load all front-space relations at weight w as a flat {dependent: {indep: coef}} dict."""
    assert (w > 0)
    with open(f'{relpath}/{prefix}', 'rt') as f:
        return {k: v for j in getFrel_eqs(f, w, prefix=prefix) for k, v in rel_to_dict(j, False).items() if k}

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
    """Load front-space dict for a given form."""
    loaded = {}
    if fForm == 'P':
        loaded["fspace"], loaded["fspace_flip"] = get_perm_fspace(fweight)
    elif fForm == 'R':
        loaded["fspace"], loaded["fspace_flip"] = get_rest_fspace(fweight)
        loaded["f_rels"] = get_frels(fweight, relpath)
    elif fForm == 'PD':
        from aiamplitudes_common_dev.sewing_matrix_tools.coproduct_utils import Pdualspace
        loaded["fspace"] = Pdualspace(fweight, 'front')
        loaded["fspace_flip"] = None
    else:
        raise ValueError
    return loaded


def _load_bspace(bForm, bweight):
    """Load back-space dict for a given form."""
    loaded = {}
    if bForm == 'P':
        loaded["bspace"], loaded["bspace_flip"] = get_perm_bspace(bweight)
    elif bForm == 'R':
        loaded["bspace"], loaded["bspace_flip"] = get_rest_bspace(bweight)
        loaded["b_rels"] = get_brels(bweight, relpath)
    elif bForm == 'PD':
        from aiamplitudes_common_dev.sewing_matrix_tools.coproduct_utils import Pdualspace
        loaded["bspace"] = Pdualspace(bweight, 'back')
        loaded["bspace_flip"] = None
    else:
        raise ValueError
    return loaded


_preload_cache = {}

def preload_fbspaces(seam, fForm, fweight, bForm, bweight):
    """Pre-load F/B-space basis dicts and relations into a single dict for batch processing.

    Results are cached for reuse across calls with the same arguments.
    """
    cache_key = (seam, fForm, fweight, bForm, bweight)
    if cache_key in _preload_cache:
        return _preload_cache[cache_key]

    loaded = dict()

    if seam == 'front':
        loaded.update(_load_fspace(fForm, fweight))
    elif seam == 'back':
        loaded.update(_load_bspace(bForm, bweight))
    elif seam == 'all':
        loaded.update(_load_fspace(fForm, fweight))
        loaded.update(_load_bspace(bForm, bweight))
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
    filepath = f'{relpath}/{filename}'
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

