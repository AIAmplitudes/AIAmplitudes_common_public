from AIAmplitudes.file_readers import relpath,convert
from AIAmplitudes.polynomial_utils import polynom_convert, get_runpolynomials, get_polynomialcoeffs
from AIAmplitudes.fbspaces import get_frels,get_brels,get_perm_fspace,get_perm_bspace
from AIAmplitudes.fbspaces import get_rest_fspace,get_rest_bspace

def Phi2Symb(L, type=None):
    if not type:
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
            print("cannot encode uncompressed!")
            raise ValueError
        if L < 7:
            symb = convert(f'{relpath}/EZ_symb_quad_new_norm', L, "quad")
        elif L == 7: symb = convert(f'{relpath}/EZ7_symb_quad_new_norm', L, "quad")
        return symb
    elif type == "oct":
        if L < 8:
            symb = convert(f'{relpath}/EZ_symb_oct_new_norm', L, "oct")
        else:
            symb = convert(f'{relpath}/EZ8_symb_oct_new_norm', L, "oct")
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

def brels():
    get_brels(relpath)

def frels():
    get_frels(relpath)

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