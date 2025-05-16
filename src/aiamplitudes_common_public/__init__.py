from aiamplitudes_common_public.download_data import relpath
from aiamplitudes_common_public.file_readers import convert,get_relpermdict
from aiamplitudes_common_public.polynomial_utils import polynom_convert, get_runpolynomials, get_polynomialcoeffs
from aiamplitudes_common_public.fbspaces import get_frels,get_brels,get_perm_fspace,get_perm_bspace
from aiamplitudes_common_public.fbspaces import get_rest_fspace,get_rest_bspace
from aiamplitudes_common_public.rels_utils import alphabet,quad_prefix


# fixed alphabet
def Phi2Symb(L, type=None):
    if not type or type == "full":
        if L > 6:
            print("cannot encode uncompressed!")
            raise ValueError
        if L==6:
            symb = convert(f'{relpath}/EZ6_symb_new_norm', L)
        else:
            symb  = convert(f'{relpath}/EZ_symb_new_norm',L)
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
    if L==6:
        symb = convert(f'{relpath}/EE33_6_symb_new_norm', L)
    else:
        symb  = convert(f'{relpath}/EE33_symb_new_norm',L)
    return symb

def Phi2Symbs():
    return {L:Phi2Symb(L) for L in [1,2,3,4,5,6]}

def Phi3Symbs():
    return {L:Phi3Symb(L) for L in [1,2,3,4,5,6]}

def runpolynomials(type=None):
    if "coeffs" in type:
        return get_polynomialcoeffs(type)
    else:
        return get_runpolynomials()

def br_rels(w,mydir=relpath):
    return get_brels(w,mydir)

def fr_rels(w,mydir=relpath):
    return get_frels(w,mydir)

def fp_1l_rels(w,mydir=relpath):
    return get_relpermdict(mydir, w, "front", "oneletter")

def fp_2l_rels(w,mydir=relpath):
    return get_relpermdict(mydir, w, "front", "twoletter")

def bp_1l_rels(w,mydir=relpath):
    return get_relpermdict(mydir, w, "back", "oneletter")

def bp_2l_rels(w,mydir=relpath):
    return get_relpermdict(mydir, w, "back", "twoletter")

def fspace(w,rp="p"):
    if rp == "p": return get_perm_fspace(w)[0]
    elif rp == "r": return get_rest_fspace(w)[0]
    else: return

def bspace(w,rp="p"):
    if rp == "p": return get_perm_bspace(w)[0]
    elif rp == "r": return get_rest_bspace(w)[0]
    else: return

def fspace_flip(w,rp="p"):
    if rp == "p": return get_perm_fspace(w)[1]
    elif rp == "r": return get_rest_fspace(w)[1]
    else: return

def bspace_flip(w,rp="p"):
    if rp == "p": return get_perm_bspace(w)[1]
    elif rp == "r": return get_rest_bspace(w)[1]
    else: return


