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
from aiamplitudes_common_public.fbspaces import get_frels,get_brels,get_perm_fspace,get_perm_bspace
from aiamplitudes_common_public.uncompressor import expand_symb
from aiamplitudes_common_public.fbspaces import get_rest_fspace,get_rest_bspace, get_perm_fspace_wt6
from aiamplitudes_common_public.rels_utils import alphabet,quad_prefix
from aiamplitudes_common_public.uncompressor import UnQuad, UnQuadLoop, UnQuadTerm
from aiamplitudes_common_public.uncompressor import UnOct, UnOctLoop, UnOctTerm
from aiamplitudes_common_public.lazy_symbol import LazySymbol

phi2_backspace_file = 'phi2multifinal_E'
phi3_backspace_file = 'phi3multifinal_E'
frontspace_file = 'multiinitial_E'

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

def br_rels(w, optrace='phi2', mydir=relpath):
    """Load back-space (restrictive) relations at weight w."""
    return get_brels(w, mydir, prefix=_backspace_files[optrace])

def fr_rels(w,mydir=relpath):
    """Load front-space (restrictive) relations at weight w."""
    return get_frels(w,mydir, prefix=frontspace_file)

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

def fspace(w,rp="P"):
    """Load front-space basis dict at weight w. rp='P' for permissive, 'R' for restrictive, 'PD' for dual permissive."""
    if rp == "PD":
        from aiamplitudes_common_dev.sewing_matrix_tools.coproduct_utils import Pdualspace
        return Pdualspace(w, 'front')
    if w == 6 and rp == "P":
        return get_perm_fspace_wt6()
    else:
        if w>0:
            if rp == "P": return get_perm_fspace(w)[0]
            elif rp == "R": return get_rest_fspace(w)[0]
            else: return
        else: return

def bspace(w, rp="P", optrace='phi2'):
    """Load back-space basis dict at weight w. rp='P' for permissive, 'R' for restrictive, 'PD' for dual permissive."""
    if rp == "PD":
        from aiamplitudes_common_dev.sewing_matrix_tools.coproduct_utils import Pdualspace
        return Pdualspace(w, 'back', optrace=optrace)
    if w>0:
        if rp == "P": return get_perm_bspace(w)[0]
        elif rp == "R": return get_rest_bspace(w, prefix=_backspace_files[optrace])[0]
        else: return
    else: return

def fspace_flip(w,rp="P"):
    """Load front-space reverse-lookup dict (letter string -> basis name) at weight w."""
    if rp == "P": return get_perm_fspace(w)[1]
    elif rp == "R": return get_rest_fspace(w)[1]
    else: return

def bspace_flip(w, rp="P", optrace='phi2'):
    """Load back-space reverse-lookup dict (letter string -> basis name) at weight w."""
    if rp == "P": return get_perm_bspace(w)[1]
    elif rp == "R": return get_rest_bspace(w, prefix=_backspace_files[optrace])[1]
    else: return


