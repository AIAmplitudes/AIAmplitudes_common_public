"""
Public API for the AIAmplitudes project.

Provides convenience functions for loading Phi2/Phi3 scattering amplitude symbols,
accessing front/back space bases and relations, decompressing quad/oct data,
and working with ray polynomials.

The 6-letter alphabet {a,b,c,d,e,f} encodes the kinematics of the 3-point
form factor. Symbols are dicts mapping letter-string keys to integer coefficients.
"""

import os
from numbers import Integral

from aiamplitudes_common_public.download_data import relpath
from aiamplitudes_common_public.file_readers import convert,get_relpermdict
from aiamplitudes_common_public.polynomial_utils import polynom_convert, get_runpolynomials, get_polynomialcoeffs
from aiamplitudes_common_public.fbspaces import get_frels,get_brels,get_perm_fspace,get_perm_bspace
from aiamplitudes_common_public.uncompressor import expand_symb
from aiamplitudes_common_public.fbspaces import get_rest_fspace,get_rest_bspace, get_perm_fspace_wt6
from aiamplitudes_common_public.fbspaces import (
    _basis_cache, basis_sizes, clear_basis_cache, _basis_filepath,
    _relation_count, relation_counts,
)
from aiamplitudes_common_public.rels_utils import alphabet,quad_prefix
from aiamplitudes_common_public.uncompressor import UnQuad, UnQuadLoop, UnQuadTerm
from aiamplitudes_common_public.uncompressor import UnOct, UnOctLoop, UnOctTerm
from aiamplitudes_common_public.lazy_symbol import LazySymbol

phi2_backspace_file = 'phi2multifinal_E'
phi3_backspace_file = 'phi3multifinal_E'
frontspace_file = 'multiinitial_E'
frontspace_perm_file = 'frontspace'
backspace_perm_file = 'backspace'

_backspace_files = {
    'phi2': phi2_backspace_file,
    'phi3': phi3_backspace_file,
}


def Phi2Symb(L, type=None, uncompress = True):
    """Load the Phi2 (tr(F^2)) form factor symbol at loop order L.

    Args:
        L: Loop order (1-8). Word length = 2L for full format;
           quad/oct keys have shorter stems plus compressed tokens.
        type: None/'full' for expanded, 'quad' for 4-letter compressed,
              'oct' for 8-letter compressed.
        uncompress: If True and type is None, auto-expands quad/oct data.
                    L=7 auto-expands from quad (slow).
                    L=8 requires explicit type='oct'; type=None raises ValueError.
    Returns:
        Dict mapping letter-string keys to integer coefficients.
    """
    if not uncompress and not type:
        if L == 7: type="quad"
        elif L== 8: type = "oct"
        else: type = "full"

    if not type or type == "full":
        if L< 6:
            symb  = convert(f'{relpath}/EZ_symb_new_norm',L)
        elif L==6:
            symb = convert(f'{relpath}/EZ6_symb_new_norm', L)
        elif L == 7 and uncompress:
            print("Expanding! This will take a while and is not recommended!")
            symb = expand_symb(convert(f'{relpath}/EZ7_symb_quad_new_norm', L, "quad"))
        else:
            raise ValueError
        return symb
    elif type == "quad":
        if L < 2:
            print("cannot encode quad!")
            raise ValueError
        if L < 7:
            symb = convert(f'{relpath}/EZ_symb_quad_new_norm', L, "quad")
        elif L == 7: symb = convert(f'{relpath}/EZ7_symb_quad_new_norm', L, "quad")
        else: raise ValueError
        return symb
    elif type == "oct":
        if L < 4:
            print("cannot encode oct!")
            raise ValueError
        if L < 8:
            symb = convert(f'{relpath}/EZ_symb_oct_new_norm', L, "oct")
        elif L==8:
            symb = convert(f'{relpath}/EZ8_symb_oct_new_norm', L, "oct")
        else: raise ValueError
        return symb
    else: return

def Phi3Symb(L):
    """Load the Phi3 (tr(F^3)) form factor symbol at loop order L (1-6)."""
    if L==6:
        symb = convert(f'{relpath}/EE33_6_symb_new_norm', L)
    else:
        symb  = convert(f'{relpath}/EE33_symb_new_norm',L)
    return symb

def Phi2Symbs():
    """Load Phi2 symbols for all available loop orders (1-6) as a dict keyed by L."""
    return {L:Phi2Symb(L) for L in [1,2,3,4,5,6]}

def Phi3Symbs():
    """Load Phi3 symbols for all available loop orders (1-6) as a dict keyed by L."""
    return {L:Phi3Symb(L) for L in [1,2,3,4,5,6]}

def runpolynomials(type=None):
    """Load d-ray polynomials. type is required: pass 'coeffs' or 'coeffs_enc' for
    coefficient form, or any other string for the raw polynomial dict."""
    if "coeffs" in type:
        return get_polynomialcoeffs(type)
    else:
        return get_runpolynomials()

def br_rels(w, optrace='phi2', mydir=relpath, *, prefix=None):
    """Load cached back relations. Treat the returned dict as read-only."""
    if prefix is None:
        prefix = _backspace_files[optrace]
    return get_brels(w, mydir, prefix=prefix)

def _br_referenced_indeps(w, optrace='phi2', mydir=relpath):
    """Return the letter words used on the right-hand side of ``br_rels``."""
    return {
        word
        for expansion in br_rels(w, optrace, mydir).values()
        for word in expansion
        if word is not None
    }

def br_mand(w, optrace='phi2', mydir=relpath):
    """Return mandatory restrictive back-space elements at weight w.

    The result is a ``{BR_name: letter_word}`` dict. Mandatory elements are
    basis representatives that do not occur on the right-hand side of any
    back-space relation.
    """
    referenced = _br_referenced_indeps(w, optrace, mydir)
    return {
        name: word
        for name, word in bspace(w, "R", optrace=optrace).items()
        if word not in referenced
    }

def br_indep(w, optrace='phi2', mydir=relpath):
    """Return referenced independent restrictive back-space elements at weight w.

    The result is a ``{BR_name: letter_word}`` dict containing the basis
    representatives that occur on the right-hand side of ``br_rels(w)``.
    """
    referenced = _br_referenced_indeps(w, optrace, mydir)
    return {
        name: word
        for name, word in bspace(w, "R", optrace=optrace).items()
        if word in referenced
    }

def br_dep(w, optrace='phi2', mydir=relpath):
    """Return dependent back-space elements and their expansions at weight w.

    This has the same ``{dependent_word: {independent_word: coefficient}}``
    representation as ``br_rels(w)`` and includes relations that vanish.
    """
    return br_rels(w, optrace, mydir)

def fr_rels(w, mydir=relpath, *, prefix=None):
    """Load cached front relations. Treat the returned dict as read-only."""
    if prefix is None:
        prefix = frontspace_file
    return get_frels(w, mydir, prefix=prefix)


def relation_count(w, space='front', *, optrace='phi2', mydir=relpath, prefix=None):
    """Return the cached count of restrictive relations in the selected data.

    Counts parsed dependent-word entries, including zero relations, rather
    than matrix rank. The relations are loaded and counted only once per
    source file, space, and weight. Weight zero has no relations. An explicit
    prefix inspects another file without changing the default selection.
    """
    if space not in {'front', 'back'}:
        raise ValueError("space must be 'front' or 'back'")
    if space == 'back' and optrace not in _backspace_files:
        raise ValueError("optrace must be 'phi2' or 'phi3'")
    if prefix is None:
        prefix = frontspace_file if space == 'front' else _backspace_files[optrace]
    return _relation_count(w, space, mydir, prefix)

def _fr_referenced_indeps(w, mydir=relpath):
    """Return the letter words used on the right-hand side of ``fr_rels``."""
    return {
        word
        for expansion in fr_rels(w, mydir).values()
        for word in expansion
        if word is not None
    }

def fr_mand(w, mydir=relpath):
    """Return mandatory restrictive front-space elements at weight w.

    The result is a ``{FR_name: letter_word}`` dict.  Mandatory elements are
    basis representatives that do not occur on the right-hand side of any
    front-space relation.
    """
    referenced = _fr_referenced_indeps(w, mydir)
    return {name: word for name, word in fspace(w, "R").items() if word not in referenced}

def fr_indep(w, mydir=relpath):
    """Return referenced independent restrictive front-space elements at weight w.

    The result is a ``{FR_name: letter_word}`` dict containing the basis
    representatives that occur on the right-hand side of ``fr_rels(w)``.
    """
    referenced = _fr_referenced_indeps(w, mydir)
    return {name: word for name, word in fspace(w, "R").items() if word in referenced}

def fr_dep(w, mydir=relpath):
    """Return dependent front-space elements and their expansions at weight w.

    This has the same ``{dependent_word: {independent_word: coefficient}}``
    representation as ``fr_rels(w)`` and includes relations that vanish.
    """
    return fr_rels(w, mydir)

def fp_1l_rels(w,mydir=relpath):
    """Load front-space permissive 1-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "front", "oneletter")

def fp_2l_rels(w,mydir=relpath):
    """Load front-space permissive 2-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "front", "twoletter")

def bp_1l_rels(w, mydir=relpath):
    """Load back-space permissive 1-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "back", "oneletter")

def bp_2l_rels(w, mydir=relpath):
    """Load back-space permissive 2-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "back", "twoletter")

def br_zeros(w, optrace='phi2', mydir=relpath):
    """Return list of back-space elements that vanish at weight w."""
    return [k for k, v in get_brels(w, mydir, prefix=_backspace_files[optrace]).items() if v == {None: 0}]

def fr_zeros(w,mydir=relpath):
    """Return list of front-space elements that vanish at weight w."""
    return [k for k, v in get_frels(w, mydir, prefix=frontspace_file).items() if v == {None: 0}]

def br_nzrels(w, optrace='phi2', mydir=relpath):
    """Return non-zero back-space relations at weight w (dependent -> {indep: coef})."""
    return {k: v for k, v in get_brels(w, mydir, prefix=_backspace_files[optrace]).items() if v != {None: 0}}

def fr_nzrels(w,mydir=relpath):
    """Return non-zero front-space relations at weight w (dependent -> {indep: coef})."""
    return {k: v for k, v in get_frels(w, mydir, prefix=frontspace_file).items() if v != {None: 0}}

def configure_basis(*, front_rest=None, back_rest=None, front_perm=None,
                    back_perm=None, optrace='phi2'):
    """Select default basis files for subsequent loads and dimension lookups.

    Paths can be absolute or relative to the data directory. Restrictive files
    also supply the associated relations. Unspecified files stay unchanged.
    All basis-dependent caches are cleared. Existing dictionaries held by the
    caller remain snapshots of the previously loaded basis.
    """
    global frontspace_file, frontspace_perm_file, backspace_perm_file
    if optrace not in _backspace_files:
        raise ValueError("optrace must be 'phi2' or 'phi3'")
    paths = [None if value is None else os.fspath(value)
             for value in (front_rest, back_rest, front_perm, back_perm)]
    front_rest, back_rest, front_perm, back_perm = paths
    if front_rest is not None:
        frontspace_file = front_rest
    if back_rest is not None:
        _backspace_files[optrace] = back_rest
    if front_perm is not None:
        frontspace_perm_file = front_perm
    if back_perm is not None:
        backspace_perm_file = back_perm
    clear_basis_cache()


def _basis_key(w, space, rp, optrace, prefix):
    if not isinstance(w, Integral) or isinstance(w, bool) or w < 0:
        raise ValueError("weight must be a nonnegative integer")
    if space not in {'front', 'back'}:
        raise ValueError("space must be 'front' or 'back'")
    if rp not in {'R', 'P', 'PD'}:
        raise ValueError("rp must be 'R', 'P', or 'PD'")
    if space == 'back' and optrace not in _backspace_files:
        raise ValueError("optrace must be 'phi2' or 'phi3'")
    if rp == 'PD' and prefix is not None:
        raise ValueError("Select the restrictive source with configure_basis for PD")
    special_front = (space == 'front' and rp == 'P' and w == 6
                     and prefix is None and frontspace_perm_file == 'frontspace')
    if prefix is None:
        if space == 'front':
            prefix = (frontspace_file if rp in {'R', 'PD'} else
                      'wt6_242_indep_symbols' if special_front else frontspace_perm_file)
        else:
            prefix = _backspace_files[optrace] if rp in {'R', 'PD'} else backspace_perm_file
    source = _basis_filepath(prefix)
    key = (space, rp, optrace if space == 'back' else None, source, int(w), special_front)
    return key


def _load_basis(key):
    if key in _basis_cache:
        return _basis_cache[key]
    space, rp, optrace, source, w, special_front = key
    tag = ('F' if space == 'front' else 'B') + rp
    if w == 0:
        basis = {f'{tag}_0_1': '' if rp == 'R' else {'': 1}}
    elif rp == 'PD':
        from aiamplitudes_common_dev.sewing_matrix_tools.coproduct_utils import Pdualspace
        basis = Pdualspace(w, space, optrace=optrace or 'phi2')
    elif special_front:
        basis = get_perm_fspace_wt6(filename=source)
    else:
        loader = {
            ('front', 'P'): get_perm_fspace,
            ('front', 'R'): get_rest_fspace,
            ('back', 'P'): get_perm_bspace,
            ('back', 'R'): get_rest_bspace,
        }[space, rp]
        basis, flip = loader(w, prefix=source)
        _basis_cache[key] = (basis, flip)
        basis_sizes[key] = len(basis)
        return _basis_cache[key]
    if rp == 'R':
        flip = {word: name for name, word in basis.items()}
    elif rp == 'PD':
        flip = None
    else:
        flip = {}
        for name, expansion in basis.items():
            for word, coefficient in expansion.items():
                flip.setdefault(word, {})[name] = coefficient
    _basis_cache[key] = (basis, flip)
    basis_sizes[key] = len(basis)
    return _basis_cache[key]


def basis_size(w, space='front', rp='R', *, optrace='phi2', prefix=None):
    """Return the cached size of a basis, loading it once if necessary.

    The cache distinguishes weight, front/back, R/P/PD, operator trace, and
    source file. No file reads or counting occur on subsequent size lookups.
    """
    key = _basis_key(w, space, rp, optrace, prefix)
    if key not in basis_sizes:
        _load_basis(key)
    return basis_sizes[key]


def fspace(w, rp="P", *, prefix=None):
    """Load and cache the selected front basis. Treat the returned dict as read-only.

    An optional prefix loads a different file without changing the defaults;
    use configure_basis(front_rest=... / front_perm=...) to select defaults.
    """
    return _load_basis(_basis_key(w, 'front', rp, 'phi2', prefix))[0]


def bspace(w, rp="P", optrace='phi2', *, prefix=None):
    """Load and cache the selected back basis. Treat the returned dict as read-only."""
    return _load_basis(_basis_key(w, 'back', rp, optrace, prefix))[0]


def fspace_flip(w, rp="P", *, prefix=None):
    """Return the cached reverse lookup for the selected front basis."""
    return _load_basis(_basis_key(w, 'front', rp, 'phi2', prefix))[1]


def bspace_flip(w, rp="P", optrace='phi2', *, prefix=None):
    """Return the cached reverse lookup for the selected back basis."""
    return _load_basis(_basis_key(w, 'back', rp, optrace, prefix))[1]
