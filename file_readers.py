#utility functions block

def fast_convert(loop, filename):
    dev = re.split(":=|SB\(|\)", re.sub('[,*]', '', readESymb_simple(loop, filename)))[1:-1]
    keys = dev[1::2]
    values = [int(re.sub('[+-]$', t[0] + '1', t)) for t in dev[0::2]]
    out_dict = {}
    for k, v in zip(keys, values):
        out_dict[k] = v
    return out_dict

def convert(loop, filename, quad, ae, aef, octuples=False, orbit=""):
    if octuples:
        base = readESymb(loop, filename, False, ae, aef, True)[:-2]
        base = re.sub(' ', '', base)
        print(base[-10:])
        t = re.split(":=\[|\),|\)\]", base)[1:]
        if len(t[-1]) == 0: t = t[:-1]
        print(len(t), "elements")
        s = [re.split(":=|SB\(|\)", re.sub('[, *]', '', tt)) for tt in t]
        prefix = [f"v{i:02d}" for i in range(93)]
        dev = []
        for i, ss in enumerate(s):
            for j, tt in enumerate(ss[1::2]):
                s[i][1 + 2 * j] = prefix[i] + tt
            dev += s[i]
        print(len(dev), "terms")
        if len(dev[-1]) == 0: dev = dev[:-1]
        print(dev[-10:])
    elif quad:
        base = readESymb(loop, filename, True, ae, aef)[:-2]
        base = re.sub(' ', '', base)
        print(base[-10:])
        t = re.split(":=\[|\),|\)\]", base)[1:]
        if len(t[-1]) == 0: t = t[:-1]
        print(len(t), "elements")
        s = [re.split(":=|SB\(|\)", re.sub('[, *]', '', tt)) for tt in t]
        prefix = list("abcdefgh")
        dev = []
        for i, ss in enumerate(s):
            for j, tt in enumerate(ss[1::2]):
                s[i][1 + 2 * j] = prefix[i] + tt
            dev += s[i]
        print(len(dev), "terms")
        if len(dev[-1]) == 0: dev = dev[:-1]
        print(dev[-10:])
    else:
        dev = re.split(":=|SB\(|\)", re.sub('[,*]', '', readESymb(loop, filename, False, ae, aef)))[1:-1]
    keys = dev[1::2]
    values = [int(re.sub('[+-]$', t[0] + '1', t)) for t in dev[0::2]]
    out_dict = {}
    for k, v in zip(keys, values):
        if orbit == "single_quad":
            # only get one quad
            for elem in quad_list:
                if k.endswith(elem):
                    out_dict[k] = v
        elif orbit == "single_triple":
            # only get one triple
            for elem in triple_list:
                if k.startswith(elem):
                    out_dict[k] = v
        elif orbit == "single_quadtrip":
            # only get one quad and one triple
            for elem in triple_list:
                for elem2 in quad_list:
                    if k.startswith(elem) and k.endswith(elem2):
                        out_dict[k] = v
        else:
            out_dict[k] = v
    return out_dict

def readESymb(loop, name, quad, ae, aef, octuples=False):
    assert os.path.isfile(name)
    res = ''
    if octuples:
        prefix = 'Esymboct'
    elif quad:
        prefix = 'Esymbquad'
    elif ae:
        prefix = 'Eae'
    elif aef:
        prefix = 'Eaef'
    else:
        prefix = 'Esymb'

    with open(name, 'rt') as f:
        reading_form = False
        for line in f:
            if not reading_form:
                if not line.startswith(prefix + '[' + str(loop) + ']'): continue
                res = ''
                reading_form = True
            res += line[:-2] if line[-2] == '\\' else line[:-1]
            if line[-2] in [":", ";"]:
                break
    return res
    
def readESymb_simple(loop, name):
    assert os.path.isfile(name)
    res = ''
    prefix = 'Esymb'

    with open(name, 'rt') as f:
        reading_form = False
        for line in f:
            if not reading_form:
                if not line.startswith(prefix + '[' + str(loop) + ']'): continue
                res = ''
                reading_form = True
            res += line[:-2] if line[-2] == '\\' else line[:-1]
            if line[-2] in [":", ";"]:
                break
    return res


