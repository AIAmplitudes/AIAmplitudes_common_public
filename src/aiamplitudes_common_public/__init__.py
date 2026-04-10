"""
Public API for the AIAmplitudes project.

Provides convenience functions for loading Phi2/Phi3 scattering amplitude symbols,
accessing front/back space bases and relations, decompressing quad/oct data,
and working with ray polynomials.

The 6-letter alphabet {a,b,c,d,e,f} encodes the kinematics of the 3-point
form factor. Symbols are dicts mapping letter-string keys to integer coefficients.
"""

from aiamplitudes_common_public.download_data import relpath
from aiamplitudes_common_public.file_readers import convert,get_relpermdict
from aiamplitudes_common_public.polynomial_utils import polynom_convert, get_runpolynomials, get_polynomialcoeffs
from aiamplitudes_common_public.fbspaces import get_frels,get_brels,get_perm_fspace,get_perm_bspace, expand_symb
from aiamplitudes_common_public.fbspaces import get_rest_fspace,get_rest_bspace
from aiamplitudes_common_public.rels_utils import alphabet,quad_prefix
from aiamplitudes_common_public.uncompressor import UnQuad, UnQuadLoop, UnQuadTerm
from aiamplitudes_common_public.uncompressor import UnOct, UnOctLoop, UnOctTerm
from aiamplitudes_common_public.lazy_symbol import LazySymbol


def Phi2Symb(L, type=None, uncompress = True):
    """Load the Phi2 (tr(F^2)) form factor symbol at loop order L.

    Args:
        L: Loop order (1-8). Word length = 2L.
        type: None/'full' for expanded, 'quad' for 4-letter compressed,
              'oct' for 8-letter compressed.
        uncompress: If True and type is None, auto-expands quad/oct data.
                    L=7 auto-expands from quad (slow). L=8 has no auto-expand.
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
    """Load d-ray polynomials. Pass type='coeffs' or 'coeffs_enc' for coefficient form."""
    if "coeffs" in type:
        return get_polynomialcoeffs(type)
    else:
        return get_runpolynomials()

def br_rels(w,mydir=relpath):
    """Load back-space (restrictive) relations at weight w."""
    return get_brels(w,mydir)

def fr_rels(w,mydir=relpath):
    """Load front-space (restrictive) relations at weight w."""
    return get_frels(w,mydir)

def fp_1l_rels(w,mydir=relpath):
    """Load front-space permissive 1-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "front", "oneletter")

def fp_2l_rels(w,mydir=relpath):
    """Load front-space permissive 2-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "front", "twoletter")

def bp_1l_rels(w,mydir=relpath):
    """Load back-space permissive 1-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "back", "oneletter")

def bp_2l_rels(w,mydir=relpath):
    """Load back-space permissive 2-letter coproduct relations at weight w."""
    return get_relpermdict(mydir, w, "back", "twoletter")

def fspace(w,rp="P"):
    """Load front-space basis dict at weight w. rp='P' for permissive, 'R' for restrictive."""
    if w>0:
        if rp == "P": return get_perm_fspace(w)[0]
        elif rp == "R": return get_rest_fspace(w)[0]
        else: return
    else: return

def bspace(w,rp="P"):
    """Load back-space basis dict at weight w. rp='P' for permissive, 'R' for restrictive."""
    if w>0:
        if rp == "P": return get_perm_bspace(w)[0]
        elif rp == "R": return get_rest_bspace(w)[0]
        else: return
    else: return

def fspace_flip(w,rp="P"):
    """Load front-space reverse-lookup dict (letter string -> basis name) at weight w."""
    if rp == "P": return get_perm_fspace(w)[1]
    elif rp == "R": return get_rest_fspace(w)[1]
    else: return

def bspace_flip(w,rp="P"):
    """Load back-space reverse-lookup dict (letter string -> basis name) at weight w."""
    if rp == "P": return get_perm_bspace(w)[1]
    elif rp == "R": return get_rest_bspace(w)[1]
    else: return


