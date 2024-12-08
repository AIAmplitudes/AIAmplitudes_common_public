from file_readers import convert

def getSymbs(type,L):
    if type=='phi2': mystr='EZ'
    elif type=='phi3':mystr='EE33_'
    else:
        print("Error! Unknown type!")
        raise ValueError
    if L==6:
        symb = convert(L,f'../data/{mystr}6_symb_new_norm')
    else:
        newstr='EZ_'  if type=='phi2' else mystr
        symb  = convert(L,f'../data/{newstr}symb_new_norm')

    return symb

symbols={L:getSymbs('phi2', L) for L in [1,2,3,4,5,6]}
symbols_phi3={L:getSymbs('phi3', L) for L in [1,2,3,4,5,6]}