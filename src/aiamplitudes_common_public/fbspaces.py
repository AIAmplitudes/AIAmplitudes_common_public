import re
import os
from fractions import Fraction
from aiamplitudes_common_public.file_readers import readSymb, readFile, SB_to_dict
from aiamplitudes_common_public.download_data import relpath

B_number= [0, 3, 6, 12, 24, 45, 85, 155, 279, None ] #<- dim_back
F_number= [0, 3, 9, 21, 48, 108, 246, 555, 1251, None ]

NB_rels= [0, 3, 12, 24, 48, 99, 185, 355, 651, None ] #<-num_bspace_rels
NF_rels= [0, 3, 9, 33] #<-num_fspace_rels

bspacenames = {1: 'singleindep3',
               2: 'doubleindep6',
               3: 'tripleindep12',
               4: 'quadindep24',
               5: 'quintindep45',
               6: 'hexindep85',
               7: 'heptindep155',
               8: 'octindep279'}

brelnames = {1: 'singlerels3',
             2: 'doublerels12',
             3: 'triplerels24',
             4: 'quadrels48',
             5: 'quintrels99',
             6: 'hexrels185',
             7: 'heptrels355',
             8: 'octrels651'}

fspacenames = {1: 'isingleindep3',
               2: 'idoubleindep9',
               3: 'itripleindep21'}

frelnames = {1: 'isinglerels3',
             2: 'idoublerels9',
             3: 'itriplerels33'}


def get_perm_fspace(w):
    prefix='frontspace'
    assert os.path.isfile(f'{relpath}/{prefix}')
    mystr = ''.join(str.split(readSymb(f'{relpath}/{prefix}','frontspace',w)))
    newstr = re.split(":=|\[|\]", mystr)[4]
    dev = [elem + ")" if elem[-1] != ")" else elem for elem in newstr.split("),") if elem]
    basedict = {f'Fp_{w}_{i}': SB_to_dict(el) for i, el in enumerate(dev)}
    flipdict = {}
    for elem, elemdict in basedict.items():
        for term, coef in elemdict.items():
            if term not in flipdict: flipdict[term] = {}
            flipdict[term][elem] = basedict[elem][term]
    return basedict, flipdict

def get_perm_bspace(w):
    prefix = 'backspace'
    assert os.path.isfile(f'{relpath}/{prefix}')
    mystr = ''.join(str.split(readSymb(f'{relpath}/{prefix}', 'backspace', w)))
    newstr = re.split(":=|\[|\]", mystr)[4]
    dev = [elem + ")" if elem[-1] != ")" else elem for elem in newstr.split("),") if elem]
    basedict = {f'Bp_{w}_{i}': SB_to_dict(el) for i, el in enumerate(dev)}
    flipdict = {}
    for elem, elemdict in basedict.items():
        for term, coef in elemdict.items():
            if term not in flipdict: flipdict[term] = {}
            flipdict[term][elem] = basedict[elem][term]
    return basedict, flipdict


def get_rest_bspace(w):
    prefix = 'multifinal_new_norm'
    assert os.path.isfile(f'{relpath}/{prefix}')
    res=readSymb(f'{relpath}/{prefix}',str(bspacenames[w]))
    myset = {elem for elem in re.split(":=\[|E\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem}
    myd = {elem: f'Br_{w}_{i}' for i, elem in enumerate(myset)}
    flip = {f'Br_{w}_{i}': elem for i, elem in enumerate(myset)}
    return flip, myd


def get_rest_fspace(w):
    prefix='ClipFrontTriple'
    assert os.path.isfile(f'{relpath}/{prefix}')
    res=readSymb(f'{relpath}/{prefix}',str(fspacenames[w]))
    myset = {elem for elem in re.split(":=\[|SB\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem}
    myd = {elem: f'Fr_{w}_{i}' for i, elem in enumerate(myset)}
    flip = {f'Fr_{w}_{i}': elem for i, elem in enumerate(myset)}
    return flip, myd


def getBrel_eqs(f, w):
    res = readFile(f, str(brelnames[w]))
    out = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res))[1:]]
    return out


def getFrel_eqs(f, w):
    res = readFile(f, str(frelnames[w]))
    out = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res))[1:]]
    return out


def rel_to_dict(relstring, bspace=True):
    # read an F/Bspace rel as a nested dict. if the rel is
    # E(abc)=-2*E(def)+4*E(bcd), return {abc: {def:-2, bcd:4}}
    def expandcoef(c): return Fraction(c + '1') if (len(c) == 1 and not c.isnumeric()) else Fraction(c)

    if bspace:
        newstring = [elem for elem in re.split("=|E\(|\)", re.sub('#', '',
                                                                  re.sub('[,*]', '', relstring))) if elem]
    else:
        newstring = [elem for elem in re.split("=|SB\(|\)", re.sub('#', '',
                                                                  re.sub('[,*]', '', relstring))) if elem]

    if len(newstring) == 0: return {None: None}
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

def get_brels(w,relpath):
    assert (w > 0 and w < 10)
    with open(f'{relpath}/multifinal_new_norm', 'rt') as f:
        return {k: v for j in getBrel_eqs(f, w) for k, v in rel_to_dict(j).items() if k}

def get_frels(w,relpath):
    assert (w > 0  and w < 4)
    with open(f'{relpath}/ClipFrontTriple', 'rt') as f:
        return {k: v for j in getFrel_eqs(f, w) for k, v in rel_to_dict(j, False).items() if k}

def all_perm_bspaces(relpath):
    perm_bspace, perm_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        perm_bspace |= get_perm_bspace(relpath, i)
