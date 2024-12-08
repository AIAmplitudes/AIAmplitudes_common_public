import re,os
from fbspaces import get_rest_bspace,get_rest_fspace
#utility functions block

relpath='./data/'
quad_list = get_rest_bspace(relpath, 4)
oct_list = get_rest_bspace(relpath, 8)
triple_list = get_rest_fspace(relpath, 3)

def convert(filename, reptype=None, loop=None):
    #reptype: quad, oct, ae, aef, None
    if reptype in {"oct","quad"}:
        if reptype== "oct":
            base = readSymb(filename, 'ESymboct', loop)[:-2]
            prefix = [f'Br_8_{i}' for i in range(93)]
        elif reptype == "quad":
            base = readSymb(filename, 'ESymbquad', loop)[:-2]
            prefix = [f'Br_4_{i}' for i in range(8)]
        base = re.sub(' ', '', base)
        t = re.split(":=\[|\),|\)\]", base)[1:]
        if len(t[-1]) == 0: t = t[:-1]
        s = [re.split(":=|SB\(|\)", re.sub('[, *]', '', tt)) for tt in t]
        dev = []
        for i, ss in enumerate(s):
            for j, tt in enumerate(ss[1::2]):
                s[i][1 + 2 * j] = prefix[i] + tt
            dev += s[i]
        if len(dev[-1]) == 0: dev = dev[:-1]
    else:
        dev = re.split(":=|SB\(|\)", re.sub('[,*]', '',
                                            readSymb(filename, 'ESymb', loop)))[1:-1]

    keys = dev[1::2]
    values = [int(re.sub('[+-]$', t[0] + '1', t)) for t in dev[0::2]]
    out_dict = {k:v for k, v in zip(keys, values)}

    return out_dict

def readSymb(filename, prefix, loop=None):
    #read a symbol, given a filename and prefix
    assert os.path.isfile(filename)
    assert prefix in {'Esymb', 'Eaef', 'Eae', 'Esymbquad', 'Esymboct',
                      'frontspace', 'backspace','all7_new_common_factor', 'all7_sub_set', 'all6abc_list'}

    if prefix not in {'all7_new_common_factor', 'all7_sub_set', 'all6abc_list'}:
        mypref=prefix + '[' + str(loop) + ']'
    else: mypref=prefix

    with open(filename, 'rt') as f:
        # read a symbol from an opened file
        reading_form = False
        res = ''
        for line in f:
            if not reading_form:
                if not line.startswith(mypref): continue
                res = ''
                reading_form = True
            if 'space' in mypref and line.isspace(): break
            res += line[:-2] if line[-2] == '\\' else line[:-1]
            if line[-2] in [":", ";"]:
                break
    return res

def SB_to_dict(mystring):
    def to_coef(mystr):
        if mystr == '-':
            return -1
        elif mystr == '':
            return 1
        else:
            return int(mystr)

    m = mystring.replace('-', '+-').replace(",", "").split('+')
    m2 = [el.replace("(", "").replace(")", "").replace("*", "").split("SB") for el in m]
    sbdict = {elem[1]: to_coef(elem[0]) for elem in m2 if len(elem) > 1}
    return sbdict

def FBconvert(w, name):
    mystr=''.join(str.split(readSymb(filename=name, name=name, loop=w)))
    newstr=re.split(":=|\[|\]",mystr)[4]
    dev = [elem+")" if elem[-1]!=")" else elem for elem in newstr.split("),") if elem]
    return [SB_to_dict(el) for el in dev]