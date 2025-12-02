import re
import os
from fractions import Fraction
from collections import Counter
from aiamplitudes_common_public.commonclasses import Symb
from aiamplitudes_common_public.file_readers import readSymb, readFile, SB_to_dict
from aiamplitudes_common_public.download_data import relpath

B_number= [1, 3, 6, 12, 24, 45, 85, 155, 279, None ] #<- dim_back
F_number= [1, 3, 9, 21, 48, 108, 246, 555, 1251, None ]

NB_rels= [0, 3, 12, 24, 48, 99, 185, 355, 651, None ] #<-num_bspace_rels
NF_rels= [0, 3, 9, 33, 78, 180, 402, 921, 2079, None] #<-num_fspace_rels



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

fspacenames = {1 : 'isingleindep3',
               2 : 'idoubleindep9',
               3 : 'itripleindep21',
               4 : 'iquadindep48',
               5 : 'iquintindep108',
               6 : 'ihexindep246',
               7 : 'iheptindep555',
               8 : 'ioctindep1251'}


frelnames = {1: 'isinglerels3',
             2: 'idoublerels9',
             3: 'itriplerels33',
             4: 'iquadrels78',
             5: 'iquintrels180',
             6: 'ihexrels402',
             7: {'zero':'iheptrels_zerolist410', 'nonzero':'iheptrels_nonzerolist511'},
             8: 'ioctrels2079',
            }

def get_perm_fspace(w):
    prefix='frontspace'
    assert os.path.isfile(f'{relpath}/{prefix}')
    mystr = ''.join(str.split(readSymb(f'{relpath}/{prefix}','frontspace',w)))
    newstr = re.split(":=|\[|\]", mystr)[4]
    dev = [elem + ")" if elem[-1] != ")" else elem for elem in newstr.split("),") if elem]
    basedict = {f'FP_{w}_{i+1}': SB_to_dict(el) for i, el in enumerate(dev)}
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
    basedict = {f'BP_{w}_{i+1}': SB_to_dict(el) for i, el in enumerate(dev)}
    flipdict = {}
    for elem, elemdict in basedict.items():
        for term, coef in elemdict.items():
            if term not in flipdict: flipdict[term] = {}
            flipdict[term][elem] = basedict[elem][term]
    return basedict, flipdict

def get_rest_bspace(w):
    prefix = 'multifinal'
    assert os.path.isfile(f'{relpath}/{prefix}')
    res=readSymb(f'{relpath}/{prefix}',str(bspacenames[w]))
    myindeps = [elem for elem in re.split(":=\[|E\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem]
    myd = {elem: f'BR_{w}_{i+1}' for i, elem in enumerate(myindeps)}
    flip = {f'BR_{w}_{i+1}': elem for i, elem in enumerate(myindeps)}
    return flip, myd

def get_rest_bspace_OLD(w):
    prefix = 'multifinal_new_norm'
    assert os.path.isfile(f'{relpath}/{prefix}')
    res=readSymb(f'{relpath}/{prefix}',str(bspacenames[w])+' ')
    myindeps = [elem for elem in re.split(":=\[|E\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem]
    myd = {elem: f'BR_{w}_{i+1}' for i, elem in enumerate(myindeps)}
    flip = {f'BR_{w}_{i+1}': elem for i, elem in enumerate(myindeps)}
    return flip, myd

def get_rest_fspace(w):
    prefix='multiinitial'
    assert os.path.isfile(f'{relpath}/{prefix}')
    res=readSymb(f'{relpath}/{prefix}',str(fspacenames[w]))
    myindeps = [elem for elem in re.split(":=\[|SB\(|\)|\]:", re.sub('[, *]', '', res))[1:] if elem]
    myd = {elem: f'FR_{w}_{i+1}' for i, elem in enumerate(myindeps)}
    flip = {f'FR_{w}_{i+1}': elem for i, elem in enumerate(myindeps)}
    return flip, myd

def getBrel_eqs(f, w):
    res = readFile(f, str(brelnames[w]))
    out = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res))[1:]]
    return out

def getFrel_eqs(f, w):
    if w < 7:
        res = readFile(f, str(frelnames[w]))
        out = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res))[1:]]
    if w == 7:
        zeros=str(frelnames[w]['zero'])
        nonzeros=str(frelnames[w]['nonzero'])

        res_zero = readFile(f, zeros)
        res_nonzero = readFile(f, nonzeros)
        out_zero = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res_zero))[1:]]
        out_nonzero = [re.sub('\s+', '', elem) for elem in re.split(":= \[|,|\] :",
                                                        re.sub(',\s*(?=[^()]*\))', '', res_nonzero))[1:]]
        out=out_zero+out_nonzero
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
    with open(f'{relpath}/multifinal', 'rt') as f:
        return {k: v for j in getBrel_eqs(f, w) for k, v in rel_to_dict(j).items() if k}

def get_brels_OLD(w,relpath):
    assert (w > 0 and w < 10)
    with open(f'{relpath}/multifinal_new_norm', 'rt') as f:
        return {k: v for j in getBrel_eqs(f, w) for k, v in rel_to_dict(j).items() if k}

def get_frels(w,relpath):
    assert (w > 0)
    with open(f'{relpath}/multiinitial', 'rt') as f:
        return {k: v for j in getFrel_eqs(f, w) for k, v in rel_to_dict(j, False).items() if k}

def all_perm_bspaces(relpath):
    perm_bspace, perm_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        perm_bspace |= get_perm_bspace(relpath, i)

def all_rest_bspaces(relpath):
    rest_bspace, rest_bspace_keysfirst = {}, {}
    for i in range(2, 9):
        rest_bspace |= get_rest_bspace(relpath, i)


########################################################
#Functions to expand and compress symbols according to different schemes
########################################################

def expand_elem(key, opt="all", loaded=None):
    #expand. If we've already loaded the fspace and bspace dict we need, use that
    def proc_elem(myelem):
        for entry in myelem.split('@'):
            yield entry
    elem = [out for out in proc_elem(key)]
    if opt == "front":
        #fweight, bweight = int(elem[0].split('_')[1]), None
        myfspace=loaded["fspace"]
        #print(elem, myfspace[elem[0]])
        mid = ''.join(elem[1:])
        exp_dict = {f'{k}{mid}':v for k, v in myfspace[elem[0]].items()}
    elif opt == "back":
        #fweight, bweight = None, int(elem[-1].split('_')[1])
        mybspace = loaded["bspace"]
        mid=''.join(elem[:-1])
        exp_dict = {f'{mid}{k}':v for k, v in mybspace[elem[-1]].items()}
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

def preload_fbspaces(seam, fForm, fweight, bForm, bweight):
    """"
    Preload the f/b spaces and rels for a given weight combo, in order to process many terms at once
    """""
    loaded=dict()

    if seam == 'front':
        if fForm =='P':
            loaded["fspace"], loaded["fspace_flip"]= get_perm_fspace(fweight)
        elif fForm =='R':
            loaded["fspace"], loaded["fspace_flip"]= get_rest_fspace(fweight)
            loaded["f_rels"] = get_frels(fweight, relpath)
        else:
            raise ValueError
    elif seam == 'back':
        if bForm == 'P':
            loaded["bspace"], loaded["bspace_flip"]= get_perm_bspace(bweight)
        elif bForm == 'R':
            loaded["bspace"], loaded["bspace_flip"]= get_rest_bspace(bweight)
            loaded["b_rels"] = get_brels(bweight, relpath)
        else:
            raise ValueError
    elif seam == 'all':
        if fForm == 'P':
            loaded["fspace"], loaded["fspace_flip"]= get_perm_fspace(fweight)
        elif fForm == 'R':
            loaded["fspace"], loaded["fspace_flip"]= get_rest_fspace(fweight)
            loaded["f_rels"] = get_frels(fweight, relpath)
        else:
            raise ValueError

        if bForm == 'P':
            loaded["bspace"], loaded["bspace_flip"]= get_perm_bspace(bweight)
        elif bForm == 'R':
            loaded["bspace"], loaded["bspace_flip"]= get_rest_bspace(bweight)
            loaded["b_rels"] = get_brels(bweight, relpath)
        else:
            raise ValueError
    else:
        raise ValueError
    return loaded

def get_elem_compression(pattern):
    parts = pattern.split('@')

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

    if 'P' in first[0]:
        fForm = 'P'
    elif 'R' in first[0]:
        fForm = 'R'

    if 'P' in last[0]:
        bForm = 'P'
    elif 'R' in last[0]:
        bForm = 'R'

    return fForm, fweight, seamsize, bForm, bweight

def get_compression(mySymb):
    pattern=next(iter(mySymb))
    return get_elem_compression(pattern)

def expand_symb(mySymb, opt="all"):
    fForm, fweight, seamsize, bForm, bweight = get_compression(mySymb)

    if fForm ==None and bweight == None:
        print("already expanded!")
        return mySymb
    elif fweight==None:
        opt = "back"
    elif bweight==None:
        opt = "front"
    else: opt = opt

    loaded = preload_fbspaces(opt, fForm, fweight, bForm, bweight)
    res=Symb()
    for key, val in mySymb.items():
        res += expand_elem(key, opt, loaded)
    return res

def expand_rform_elem(key, opt="all", loaded=None):
    #expand. If we've already loaded the fspace and bspace dict we need, use that
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
    if seam in {"front", "all"}:
        fseam = fweight
        myfflip = loaded["fspace_flip"] if loaded else get_rest_fspace(fweight)[1]
        myf = [myfflip[''.join(key[:fseam])]]
    if seam in {"back", "all"}:
        bseam = len(key) - bweight
        mybflip = loaded["bspace_flip"] if loaded else get_rest_bspace(bweight)[1]
        myb = [mybflip[''.join(key[bseam:])]]

    if seam == "front":
        l = myf + [''.join(key[fseam:])]
    elif seam == "back":
        l = [''.join(key[:bseam])] + myb
    elif seam == "all":
        l = myf + [''.join(key[fseam:bseam])] + myb
    return '@'.join(l)

def get_related(rel, seed, slot):
    # given: a seam letter and indep (restrictive-form) bspace elem
    # i.e. c@Fr_2_3
    # AND: a nonlocal relation (integ, triple, etc.)
    # get: the other non-indep terms it is related to by this relation
    seedword = list(seed)
    return {','.join(seedword[:slot] + [i for i in elem] + seedword[slot + len(elem):]): relcoef
            for elem, relcoef in rel.items()}


def get_as_indepsum(fullword, fweight=None, bweight=None, seam=None,
                    loaded=None):
    # Given a non-indep term, get the indeps it is related to by bspace/fspace rels
    if seam == "back":
        mybflip = loaded["bspace_flip"] if loaded else get_rest_bspace(bweight)[1]
        mybrels = loaded["b_rels"] if loaded else get_brels(bweight,relpath)
        body, belem = fullword[:-bweight], ''.join(fullword[-bweight:])
        if belem in mybflip:
            reldict = {compress_rform_elem(fullword, fweight, bweight, seam, loaded): 1}
        elif belem in mybrels:
            reldict = {compress_rform_elem(body + [i for i in k], fweight, bweight, seam, loaded): v
                       for k, v in mybrels[belem].items() if k}
        else: reldict = {}

    elif seam == "front":
        myfflip = loaded["fspace_flip"] if loaded else get_rest_fspace(fweight)[1]
        myfrels = loaded["f_rels"] if loaded else get_frels(fweight,relpath)
        felem, body = ''.join(fullword[:fweight]), fullword[fweight:]
        if felem in myfflip:
            reldict = {compress_rform_elem(fullword, fweight, bweight, seam,loaded): 1}
        elif felem in myfrels:
            reldict = {compress_rform_elem([i for i in k] + body, fweight, bweight, seam, loaded): v
                       for k, v in myfrels[felem].items() if k}
        else: reldict = {}
    else:
        print("Bad seam!")
        raise ValueError
    return (reldict)












































