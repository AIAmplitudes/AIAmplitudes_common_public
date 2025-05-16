import os, re, math
from sympy import Poly, terms_gcd, gcd_list, Rational
from sympy import factorint, fraction, primerange
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication,
)
from aiamplitudes_common_public.file_readers import readSymb
from aiamplitudes_common_public.download_data import relpath

primes = list(primerange(0, 1000))

########### Encoders convert between rationals and prime factors #################
def int_to_factors(value):
    prefix = []
    if abs(value) == 1:
        factors = {1: 1}
    else:
        factors = factorint(abs(value))
    remainder = value
    for i, fac in enumerate(factors.keys()):
        if int(fac) > 0:
            remainder = int(remainder / pow(fac, factors[fac]))
            if i != 0:
                if factors[fac] != 1:
                    prefix += (f"*{fac}E{factors[fac]}")
                else:
                    prefix += (f"*{fac}")
            else:
                if factors[fac] != 1:
                    prefix += (f"{fac}E{factors[fac]}")
                else:
                    prefix += (f"{fac}")
    if remainder != 0:
        suffix = []
        w = abs(remainder)
        if w != 1:
            suffix.append(str(w))
    else:
        prefix = '0'
    prefix = ''.join(prefix)
    prefix = ('+' if value >= 0 else '-') + prefix
    return ''.join(prefix)

def frac_to_factors(value):
    val = fraction(value)
    if val[1] == '1':
        prefix = int_to_factors(val[0])
    else:
        prefix = int_to_factors(val[0]) + '/' + int_to_factors(val[1])
    prefix = ('+' if value >= 0 else '-') + prefix
    return prefix

def enc_elem(elem):
    if '/' in str(elem):
        return str(int_to_factors(elem.p)) + '/' + str(int_to_factors(elem.q))
    else:
        return str(int_to_factors(elem))

def cl(text):
    return re.sub(r'\s+', '', text)

def parse(t):
    return parse_expr(re.sub('\^','**',t),
                      transformations=standard_transformations + (implicit_multiplication,))

def polynom_convert(filename):
    input_str = readSymb(f'{relpath}/{filename}', filename)
    pattern = r'\[(.*?)\](.*?),'
    matches = re.findall(pattern, input_str)
    out_dict = {cl(re.sub(',', '', match[0])):
                    cl(re.sub('=', '', match[1])) for match in matches}
    return out_dict

def is_pow2(n):
    return (n != 0) and (n & (n - 1) == 0)

def get_runpolynomials():
    allpolys = polynom_convert('all7_new_common_factor')
    nonzeros = {k:v for k, v in allpolys.items() if v != '0'}
    unfactorable = {k:v for k, v in nonzeros.items() if '(' not in v}
    return {'all': allpolys, 'nonzero':nonzeros, 'unfactorable':unfactorable}

def get_polynomialcoeffs(type):
    mydict=get_runpolynomials()
    mycoeffs={}
    mygcds={}
    mycoeffs_nogcd={}
    mycoeffs_unenc={}
    myquadrats={}
    mylins={}
    for k,v in mydict['nonzero'].items():
        this_coeffs=Poly(parse(v).expand()).all_coeffs()
        mygcd=gcd_list(this_coeffs)
        mygcds[k]=mygcd
        if Poly(parse(v).expand()).degree('L') < 3:
            myquadrats[k] = this_coeffs
            this_coeffs = [0] + this_coeffs
            if Poly(parse(v).expand()).degree('L') < 2:
                mylins[k] = this_coeffs
                this_coeffs = [0] + this_coeffs
        myc=[c/mygcd for c in this_coeffs]
        mycoeffs_nogcd[k] = [c for c in this_coeffs]
        mycoeffs_unenc[k] = [[str(mygcd)],[c for c in myc]]
        mycoeffs[k] = [[str(mygcd)],[int_to_factors(c) for c in myc]]

    myints = {k: mycoeffs_unenc[k] for k, v in mycoeffs.items() if ('/' not in v[0][0])}
    mydivs = {k: v for k, v in mycoeffs.items() if ('/' in v[0][0])}
    # polynomials
    if type == "coeffs":
        my_all = mydivs | myints
        return {'all':my_all, 'intcoeffs':myints, 'divcoeffs':mydivs}
    elif type == "coeffs_enc":
        myints_enc = {k: [enc_elem(coef) for coef in mycoeffs_nogcd[k]] for k, v in myints.items()}
        mydivs_enc = {k: [enc_elem(coef) for coef in mycoeffs_nogcd[k]] for k, v in mydivs.items()}
        my_all_enc = mydivs_enc | myints_enc
        return {'all':my_all_enc, 'intcoeffs':myints_enc, 'divcoeffs':mydivs_enc}
    else:
        return
