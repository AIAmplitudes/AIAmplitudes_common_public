import numpy as np
from collections import Counter
from AIAmplitudes.commonclasses import Symb
import itertools
from AIAmplitudes.rels_utils import kbits,count_zeros,is_ok_phi2,is_ok_phi3,is_ok

#####################################################################################
#TODO: ADD deconcatenation & concatenation coproducts, 'cut-shuffle', other coproduct variants
#####################################################################################

def ShuffleProduct(l1, l2, cond=None):
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str)
    stufslist = [''.join(shuf) for shuf in StuffleIterator(l1, l2, r=0) if func(''.join(shuf))]
    return Symb(dict(Counter(stufslist)))

def StuffleProduct(l1, l2, cond=None):
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str)
    stufslist = []
    for r in range(min(len(l1), len(l2)) + 1):
        stufslist += [''.join(shuf) for shuf in StuffleIterator(l1, l2, r) if func(''.join(shuf))]
    return Symb(dict(Counter(stufslist)))

def UnshuffleProduct(shuffled_dict, in_dict, opt=None, cond=None, guess=False):
    # Need to be careful since L and R unshufs aren't commutative.
    # Defining a convention here: if A (LeftShuf) B = C,
    # Then B = C (LeftUnshufFirst) A, and  A = C (LeftUnshufSecond) B

    unshuf_counters = {}
    unshuf_dicts = {}
    all_possible_keys = {}
    # get a word in the shuffle output.
    for long_word in shuffled_dict:
        unshuf_counter = {}
        unshuf_dict = {}
        for k in in_dict:
            m = len(long_word)
            n = len(k)

            wc1, wc2 = long_word, k

            blank = [None] * (m - n)
            # Check that we have at least one valid way to get it with shuffles.
            # Gets: unshuffle keys and multiplicity.
            for iv in kbits(m - n, m):
                if opt == "left_first":
                    if iv[0] != 0:
                        continue
                if opt == "right_first":
                    if iv[-1] != 1:
                        continue
                if opt == "left_second":
                    if iv[0] != 1:
                        continue
                if opt == "right_second":
                    if iv[-1] != 0:
                        continue
                if opt == "leftright_first":
                    if iv[0] != 0:
                        continue
                    if iv[-1] != 1:
                        continue
                if opt == "leftright_second":
                    if iv[0] != 1:
                        continue
                    if iv[-1] != 0:
                        continue

                if count_zeros(iv) != len(k): continue
                # Fill in w1 into the iv slots
                i = 0
                is_bad = False
                unshuf_seq = []
                for j in range(len(iv)):
                    if (iv[j] == 0):
                        # if we hit a 0 and the sequences line up, we're good.
                        # otherwise, throw away this kbits and get the next one
                        if wc1[j] == wc2[i]:
                            i += 1
                        else:
                            is_bad = True
                    if (iv[j] == 1):
                        unshuf_seq += wc1[j]
                if not is_bad:
                    goodkey = ''.join(unshuf_seq)
                    if not goodkey in unshuf_counter:
                        unshuf_dict[goodkey] = k
                        unshuf_counter[goodkey] = 1
                        if not goodkey in all_possible_keys:
                            all_possible_keys[goodkey] = 1
                    else:
                        unshuf_counter[goodkey] += 1

            # for each word, store the possible unshuffles and the multiplicity
            unshuf_counters[long_word] = unshuf_counter
            unshuf_dicts[long_word] = unshuf_dict

    def get_in_coef(long_word, outkey):
        # get the coef of the input word corresponding to a given shuffle candidate,
        # times the shuffle multiplicity
        if outkey in unshuf_dicts[long_word]:
            in_coef = unshuf_counters[long_word][outkey] * in_dict[unshuf_dicts[long_word][outkey]]
        else:
            # print(unshuf_counters[long_word])
            in_coef = 0
        return in_coef

    # We know that shuffle_multiplicity_1*(in_coeff_1,shuf_coef_1)+ ... +
    # shuffle_multiplicity_n*(in_coeff_n,shuf_coef_n) = out_coeff.
    # so we need to solve a system. Construct the input matrix: [N_shuffled_words*N_unshuffle_candidates]
    M, y = [], []
    for i, long_word in enumerate(shuffled_dict):
        y.append(shuffled_dict[long_word])
        M.append([get_in_coef(long_word, key) for key in all_possible_keys.keys()])

    npM = np.array(M)
    npy = np.array(y)

    # print(npM,npy)
    outs = np.round(np.linalg.lstsq(npM, npy, rcond=None)[0]).astype(int)
    unshuf_dict = {}
    for i, key in enumerate(all_possible_keys.keys()):
        if outs[i] == 0:
            continue
        else:
            unshuf_dict[key] = outs[i]

    # Check:
    if opt == "left_first":
        mysymb = Lshufsymbs(Symb(in_dict), Symb(unshuf_dict), cond=cond)
    elif opt == "right_first":
        mysymb = Rshufsymbs(Symb(in_dict), Symb(unshuf_dict), cond=cond)
    elif opt == "left_second":
        mysymb = Lshufsymbs(Symb(unshuf_dict), Symb(in_dict), cond=cond)
    elif opt == "right_second":
        mysymb = Rshufsymbs(Symb(unshuf_dict), Symb(in_dict), cond=cond)
    elif opt == "leftright_first":
        mysymb = LRshufsymbs(Symb(in_dict), Symb(unshuf_dict), cond=cond)
    elif opt == "leftright_second":
        mysymb = LRshufsymbs(Symb(unshuf_dict), Symb(in_dict), cond=cond)
    else:
        mysymb = shufsymbs(Symb(unshuf_dict), Symb(in_dict), cond=cond)

    if mysymb != Symb(shuffled_dict):
        print("No solution!")
        if not guess: return None
    return Symb(unshuf_dict)

def LHalfShuffleProduct(l1, l2, cond=None):
    # all shuffles that begin like l1
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str)
    stufslist = [''.join(shuf) for shuf in StuffleIterator(l1, l2, r=0, opt="left") if func(''.join(shuf))]
    return Symb(dict(Counter(stufslist)))


def RHalfShuffleProduct(l1, l2, cond=None):
    # all shuffles that end like l2
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str)
    stufslist = [''.join(shuf) for shuf in StuffleIterator(l1, l2, r=0, opt="right") if func(''.join(shuf))]
    return Symb(dict(Counter(stufslist)))


def LRHalfShuffleProduct(l1, l2, cond=None):
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str)
    stufslist = [''.join(shuf) for shuf in StuffleIterator(l1, l2, r=0, opt="leftright") if func(''.join(shuf))]
    return Symb(dict(Counter(stufslist)))

def LeftPreLieShuffleProduct(l1, l2, cond=None):
    return LHalfShuffleProduct(l1, l2, cond) - RHalfShuffleProduct(l2, l1, cond)

def RightPreLieShuffleProduct(l1, l2, cond=None):
    return RHalfShuffleProduct(l1, l2, cond) - LHalfShuffleProduct(l2, l1, cond)

####################################################
#full-symbol shuffles
###################################################

def shufsymbs(symb1, symb2, cond=None):
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * ShuffleProduct(k1, k2, cond=cond)
            shufs = shufs + shuf
    return shufs


def Rshufsymbs(symb1, symb2, cond=None):
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * RHalfShuffleProduct(k1, k2, cond=cond)
            shufs = shufs + shuf
    return shufs


def Lshufsymbs(symb1, symb2, cond=None):
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * LHalfShuffleProduct(k1, k2, cond=cond)
            shufs = shufs + shuf
    return shufs


def LRshufsymbs(symb1, symb2, cond=None):
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * LRHalfShuffleProduct(k1, k2, cond=cond)
            shufs = shufs + shuf
    return shufs

def LeftPreLieShufSymbs(symb1, symb2, cond=None):
    # check
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * (LHalfShuffleProduct(k1, k2, cond=cond) - RHalfShuffleProduct(k2, k1, cond=cond))
            shufs = shufs + shuf
    return shufs

def RightPreLieShufSymbs(symb1, symb2, cond=None):
    # check
    shufs = {}
    for k1, v1 in symb1.items():
        for k2, v2 in symb2.items():
            shuf = v1 * v2 * (RHalfShuffleProduct(k1, k2, cond=cond) - LHalfShuffleProduct(k2, k1, cond=cond))
            shufs = shufs + shuf
    return shufs

##################################
#utils for the stuffle product
################################
def ShuffleProduct_overlapping_r(l1, l2, r, cond=None):
    """
    The overlapping shuffle product of the two words ``w1`` and ``w2``
    with precisely ``r`` overlaps.
    """
    if cond == "ok_phi2":
        func = is_ok_phi2
    elif cond == "ok_phi3":
        func = is_ok_phi3
    else:
        func = is_ok
    assert isinstance(l1, str) and isinstance(l2, str) and isinstance(r, int)
    _l1 = list(l1)
    _l2 = list(l2)

    return Symb(dict(Counter([''.join(shuf) for shuf in StuffleIterator(_l1, _l2, r) if func(''.join(shuf))])))

def StuffleIterator(_l1, _l2, r, opt=None):
    m = len(_l1)
    n = len(_l2)

    wc1, wc2 = _l1, _l2

    blank = [None] * (m + n - r)
    for iv in kbits(m, m + n - r):
        if opt == "left":
            if iv[0] != 1:
                continue
        if opt == "right":
            if iv[-1] != 0:
                continue
        if opt == "leftright":
            if iv[0] != 1:
                continue
            if iv[-1] != 0:
                continue
        w = blank[:]
        filled_places = []
        unfilled_places = []
        # Fill in w1 into the iv slots
        i = 0
        for j in range(len(iv)):
            if iv[j] == 1:
                w[j] = wc1[i]
                i += 1
                filled_places.append(j)
            else:
                unfilled_places.append(j)
        # Choose r of these filled places
        for subset in itertools.combinations(filled_places, r):
            places_to_fill = sorted(unfilled_places + list(subset))

            # Fill in w2 into the places
            i = 0
            res = w[:]
            for j in places_to_fill:
                if res[j] is not None:
                    res[j] = "(" + res[j] + '+' + wc2[i] + ")"
                else:
                    res[j] = wc2[i]
                i += 1
            yield [i for i in res if i]

