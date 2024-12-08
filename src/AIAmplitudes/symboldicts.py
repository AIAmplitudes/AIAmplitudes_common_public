from AIAmplitudes.file_readers import convert, relpath

def getSymbs(type,L):
    if type=='phi2': mystr='EZ'
    elif type=='phi3':mystr='EE33_'
    else:
        print("Error! Unknown type!")
        raise ValueError
    if L==6:
        symb = convert(f'{relpath}/{mystr}6_symb_new_norm', L)
    else:
        newstr='EZ_' if type=='phi2' else mystr
        symb  = convert(f'{relpath}/{newstr}symb_new_norm',L)

    return symb

symbols={L:getSymbs('phi2', L) for L in [1,2,3,4,5,6]}
phi3symbols={L:getSymbs('phi3', L) for L in [1,2,3,4,5,6]}