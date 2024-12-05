
from itertools import permutations, islice
from rels_utils import find_all

alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
dihedral_table = [list(permutations(alphabet[:3]))[i] + list(permutations(alphabet[3:]))[i] for i in
                  range(len(alphabet))]
cycle_table = [dihedral_table[i] for i in [0, 3, 4]]
flip_table = [dihedral_table[i] for i in [0, 1, 2, 5]]

##########################################################################################
#Operator Library. One-to-one operators:
###########################################################################################

def strike(key, slotset):
    # k: all slots in the operation must be within 'k' of each other
    # assume slotset is sorted, all slots are < len(key),
    # and slotsets are all within k. This is done at slot generation time.
    # can't strike slots that are not in the key!
    slotset=set(slotset)
    strike_key=''.join([k for i,k in enumerate(key) if i not in slotset])
    return strike_key

def swap(key, slotset):
    # swap letters at the given slot.
    if len(slotset) > 2: return
    swap_key=key[:slotset[0]] + key[slotset[1]] + key[(slotset[0]+1):slotset[1]] + key[slotset[0]]
    if slotset[1] != len(key): swap_key += key[(slotset[1]+1):]
    return swap_key

def unstrike(key, slotset, letterset,runlenset=None,targetloop=None):
    # insert runs of the specified letter into the key.
    # careful with lookups here!
    newkey=key
    if not runlenset:
        runlenset=[1]*len(slotset)
    for slot,let,runlen in zip(slotset,letterset,runlenset):
        newkey = newkey[:slot] + let*runlen + newkey[slot:]

    #check to make sure adding runs gets us to the correct loop
    if targetloop and(len(newkey)!=2*len(targetloop)): raise ValueError
    return newkey

def mutate(key, slotset, letterset):
    # Change the letter in the given slot to the given letter.
    keyslots=dict(zip(slotset,letterset))
    return ''.join([elem if i not in slotset else keyslots[i] for i,elem in enumerate(key)])

def subword(key,slotset):
    # get the subword that spans the two slots described,
    # INCLUDING the last letter
    if len(slotset) > 2: return
    return key[slotset[0]:slotset[1]]+key[slotset[1]]

def slot_select(key,slotset):
    return ''.join([key[i] for i in slotset])

def dihedral(key,rot_ind):
    # rot_ind is 1,2,3,4 or 5
    if rot_ind not in {1,2,3,4,5}: print(rot_ind); raise ValueError
    word_idx = [alphabet.index(l) for l in [*key]]
    image = ''.join([dihedral_table[rot_ind][idx] for idx in word_idx])
    return image

def cycle(key,rot_ind):
    #rot_ind is 3 or 4
    if rot_ind not in {3,4}: print(rot_ind); raise ValueError
    word_idx = [alphabet.index(l) for l in [*key]]
    cycle = ''.join([dihedral_table[rot_ind][idx] for idx in word_idx])
    return cycle

def flip(key,rot_ind):
    #rot_ind is 1 or 2, 5
    if rot_ind not in {1,2,5}: print(rot_ind); raise ValueError
    word_idx = [alphabet.index(l) for l in [*key]]
    flip = ''.join([flip_table[rot_ind][idx] for idx in word_idx])
    return flip

def partial_dihedral(key, slotset, rot_ind):
    if len(slotset) != 2: raise ValueError
    return key[:slotset[0]]+dihedral(key[slotset[0]:slotset[1]],rot_ind)+key[slotset[1]:]

def partial_cycle(key, slotset, rot_ind):
    if len(slotset) != 2: raise ValueError
    return key[:slotset[0]]+cycle(key[slotset[0]:slotset[1]],rot_ind)+key[slotset[1]:]

def partial_flip(key, slotset, rot_ind):
    if len(slotset) != 2: raise ValueError
    return key[:slotset[0]]+flip(key[slotset[0]:slotset[1]],rot_ind)+key[slotset[1]:]

###########################################################################
# one-to-many operators (allstrikes,etc.)
##########################################################################

def all_partial_dihedral_operator(key,slotset):
    prefix=key[:slotset[0]]
    substr=key[slotset[0]:slotset[1]]
    suffix=key[slotset[1]:]
    return [prefix+im+suffix for im in all_dihedral_operator(substr)]

def all_partial_cycle_operator(key,slotset):
    prefix=key[:slotset[0]]
    substr=key[slotset[0]:slotset[1]]
    suffix=key[slotset[1]:]
    return [prefix+im+suffix for im in all_cycle_operator(substr)]

def all_partial_flip_operator(key,slotset):
    prefix=key[:slotset[0]]
    substr=key[slotset[0]:slotset[1]]
    suffix=key[slotset[1]:]
    return [prefix+im+suffix for im in all_flip_operator(substr)]

def all_dihedral_operator(word):
    word_idx = [alphabet.index(l) for l in [*word]]
    dihedral_images = [''.join([dihedral_table[row][idx] for idx in word_idx]) for row in range(1,len(alphabet))]
    return dihedral_images

def all_cycle_operator(word):
    word_idx = [alphabet.index(l) for l in [*word]]
    cycle_images = [''.join([cycle_table[row][idx] for idx in word_idx]) for row in range(1,int(len(alphabet) / 2))]
    return cycle_images

def all_flip_operator(word):
    my_images = all_dihedral_operator(word)
    word_idx = [alphabet.index(l) for l in [*word]]
    cycle_images = [''.join([cycle_table[row][idx] for idx in word_idx]) for row in range(1,int(len(alphabet) / 2))]
    return [im for im in my_images if im not in cycle_images]

#def all_strike_operator(word):


##########################################################################
# slot selectors. Not used in instance gen, but can be in future- need to figure out
# more about how python handles passing compositions of functions
##########################################################################

def get_nth_occ_slotlookup(key, substr, n):
    #get the nth occurrence of a substring (assumes it appears >= n times)
    allapps=find_all(key, substr)
    if n > len(allapps): return None
    *_, last= islice(allapps, n)
    return last

def get_nth_occ_right_slotlookup(key, substr, n):
    #get the nth occurrence from the right of a substring (assumes it appears >= n times)
    allapps=[i for i in find_all(key, substr)]
    if n > len(allapps): return None
    last=allapps[-n]
    #*_, last= islice(allapps, n)
    return last

def get_first_nth_appearance(key, substr, n):
    if n < 1: raise ValueError
    # get the nth occurrence of a substring (assumes it appears >= n times)
    allapps = [i for i in find_all(key, substr)]
    if n > len(allapps): return None
    first = allapps[n - 1]
    return first

def get_last_nth_appearance(key, substr, n):
    if n < 1: raise ValueError
    # get the nth occurrence from the right of a substring (assumes it appears >= n times)
    allapps = [i for i in find_all(key, substr)]
    if n > len(allapps): return None
    last = allapps[-n]
    return last

def get_samerun_slots(sent, seedslot):
    thischar = sent[seedslot]
    counter = 1
    good = []
    # get all slots that are in the same "run" as the given slot
    r1, r2 = 0, 0
    # forward
    for i in range(len(sent) - seedslot + 1):
        if sent[seedslot + i] == thischar:
            continue
        else:
            r1 = seedslot + i - 1
            break
    # back
    for i in range(seedslot):
        if sent[seedslot - i] == thischar:
            continue
        else:
            r2 = seedslot - i - 1
            break
    return [r1 + i for i in range(0, r2 - r1 - 1)]


def get_runbound_slots(sent, seedslot):
    # get the slots that bookend the given "run".
    run = get_samerun_slots(sent, seedslot)
    return [run[0] - 1, run[-1] + 1]

def get_rundict(sent):
    thischar = 0
    runs = {}
    counter = 1
    for i in range(len(sent)):
        if sent[i] == thischar:
            counter += 1
        else:
            thischar = sent[i]
            if i != 0:
                if not sent[i - 1] in runs: runs[sent[i - 1]] = []
                runs[sent[i - 1]] += [(i - 1, f'r{counter}')]
            counter = 1

    if not sent[i - 1] in runs: runs[sent[i - 1]] = []
    runs[sent[i - 1]] += [(i - 1, f'r{counter}')]
    return runs

def get_first_nth_run_of(k, let, n):
    if n < 1:
        return get_rundict(k)[let][n - 1][0]

def get_last_nth_run_of(k, let, n):
    return get_rundict(k)[let][::-1][n - 1][0]

