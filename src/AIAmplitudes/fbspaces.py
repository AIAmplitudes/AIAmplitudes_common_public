import re
import os
from fractions import Fraction
from AIAmplitudes.file_readers import readSymb, SB_to_dict,relpath

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

def get_brels(relpath):
    with open(f'{relpath}/multifinal_new_norm', 'rt') as f:
        brels = {}
        for i in range(2, 9):
            brels[i] = {k: v for j in getBrel_eqs(f, i) for k, v in rel_to_dict(j).items() if k}
        return brels

def get_frels(relpath):
    with open(f'{relpath}/ClipFrontTriple', 'rt') as f:
        brels = {}
        for i in range(1, 4):
            brels[i] = {k: v for j in getFrel_eqs(f, i) for k, v in rel_to_dict(j, False).items() if k}
        return brels

def all_perm_bspaces(relpath):
    perm_bspace, perm_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        perm_bspace |= get_perm_bspace(relpath, i)

def get_rels_perm(mydir, weight, seam="front", reltype="oneletter"):
    # get the permissive rels pre-generated from seam_sewing_rels.py.
    # These are the checks generated by Lance & Garrett (lettbackrels and lettfrontrels).

    if not reltype in ["oneletter", "twoletter"]: return
    r = reltype

    if seam == "back":
        file = f"Bspace_rels_{r}"
        prefix = "sewrelsb"
    elif seam == "front":
        file = f"Fspace_rels_{r}"
        prefix = "sewrelsf"
    else:
        print("bad seam type!")
        raise ValueError


    assert os.path.isfile(f'{mydir}/{file}')
    checks = readSymb(f'{mydir}/{file}', prefix,weight)
    c = [i for i in checks.split("\'")[:-1] if i != ', ' and 'sewrels' not in i]
    return c

def readcrel(crel, w=2, seam="back"):
    def cterm_to_Fp(cterm):
        if not ',' in cterm:
            print("bad term!")
            raise ValueError
        numlet = re.split(",", cterm)
        o = []
        for i in numlet:
            if i.isnumeric():
                if seam == "front":
                    o.append(f'Fp_{w}_{i}@')
                elif seam == "back":
                    o.append(f'@Bp_{w}_{i}')
                else:
                    raise ValueError
            else:
                o.append(i)
        return ''.join(o)

    def numstr_to_num(numstr):
        if numstr == '+':
            return 1
        elif numstr == '-':
            return -1
        elif '*' in numstr:
            numstr = numstr.replace('−', '-')
            return Fraction(numstr[:-1])
            # else:return int(numstr[:-1])

    myrel = re.split("c\[|\]", crel)[:-1]
    coefs = [elem for i, elem in enumerate(myrel) if i % 2 == 0]
    terms = [elem for i, elem in enumerate(myrel) if i % 2 == 1]
    return {cterm_to_Fp(let): numstr_to_num(num) for num, let in zip(coefs, terms)}

def get_relpermdict(mydir, w, seam, reltype):
    return [readcrel(i, w, seam) for i in get_rels_perm(mydir, w, seam, reltype)]

