import itertools
import datetime

from scipy.special import binom
import random
from aiamplitudes_common_public.rels_utils import get_coeff_from_word,check_slot,find_all,alphabet,count_appearances
from aiamplitudes_common_public.commonclasses import fastRandomSampler

##########################
# generators for op_args
# note: we use tuples rather than lists for all slot combos, letter combos, etc. because tuples are hashable,
# so we can look up whether or not we've seen them before
##########################

def gen_slotsets(keylen,k_total,k_pairwise, nslots,exact=False):
    def slotcombo_generator():
        # generate all combos of n_slots that are within k of each other
        # ('within' here meaning that max(slot) - min(slot) <= k)
        # if exact, only generate those that are exactly k apart
        # if allcombos, generate all combos of
        if nslots > k_total: raise ValueError
        for i in range(0, keylen - nslots):
            kbound = min(i + k_total, keylen)
            for c in itertools.combinations(range(i, kbound), nslots):
                if (exact and (c[-1] - c[0] == k_total)) or (c[-1] - c[0] <= k_total):
                    if all([c[i+1]-c[i] <= k_pairwise for i in range(0,len(c)-1)]): yield c
    return set(slot for slot in slotcombo_generator())

def gen_kpattern_slotsets(keylen,k_pattern,exact=False):
#TODO: Finish this
    def slotcombo_generator(k_pattern):
        # generate all combos of n_slots with spacing at least (k0,k1,...)
        k_total=sum(k for k in k_pattern)
        nslots=len(k_pattern)
        if nslots > k_total: raise ValueError

        for i in range(0, keylen - nslots): #choose the start position
                                            #move a pointer from 0 to i
                                            #pop k0, update the budget.
                                            #get the next start that is within keylen-budget

            if len(k_pattern) > 1: slotcombo_generator(k_pattern[1:])  #get all combos from k_pattern[1:]
            else:
                kbound = min(i + k_total, keylen-k_total)
                for c in itertools.combinations(range(i, kbound), nslots):
                    if (exact and (c[-1] - c[0] == k_total)) or (c[-1] - c[0] <= k_total):
                        if all([c[i+1]-c[i] >= k_pattern[i] for i in range(0,len(c)-1)]):
                            yield c

    return set(slot for slot in slotcombo_generator())

def gen_lettersets(nletts):
    def lettercombo_generator():
        for c in itertools.combinations_with_replacement(alphabet, nletts):yield c
    return set(c for c in lettercombo_generator())

def gen_sumtuples(n_elems, target_sum):
    #generate all tuples of n_elems nums that sum to a target.
    #used to insert runs into a key
    def sumnum_combo_generator(n_elems, target_sum):
        if n_elems == 0: return
        # generate all sets of n_elems numbers that sum to target_sum.
        if n_elems == 1:
            yield (target_sum,)
        for i in range(1,target_sum):
            for t in sumnum_combo_generator(n_elems - 1, target_sum - i):
                yield (i,) + t
    return set(tup for tup in sumnum_combo_generator(n_elems, target_sum))

def gen_op_args(op_argdict):
    op_args=[]
    if "slots" in op_argdict:
        if "allcombos" in op_argdict["slots"]:
            op_args.append({(element,) for slotsetsize in range(2,2 * op_argdict["slots"]["loop"],2) for element in
                            gen_slotsets(2 * op_argdict["slots"]["loop"], op_argdict["slots"]["k_total"],
                                         op_argdict["slots"]["k_pairwise"], slotsetsize)})
        else:
            op_args.append(gen_slotsets(2 * op_argdict["slots"]["loop"], op_argdict["slots"]["k_total"],
                                        op_argdict["slots"]["k_pairwise"], op_argdict["slots"]["numslots"]))
    if "letters" in op_argdict:
        op_args.append(gen_lettersets(op_argdict["letters"]["numslots"]))
    if "sumtups" in op_argdict:
        op_args.append(gen_sumtuples(op_argdict["sumtups"]["numslots"],op_argdict["sumtups"]["totalmult"]))
    if "rot_ind" in op_argdict:
        op_args.append({1,2,3,4,5})
        #make tuples that hold combos of slots, letters, and other args required by the operator
    if "runs" in op_argdict:
        if "allcombos" in op_argdict["runs"]:
            #all possible sets of runs with total length of runs in set < totalmult
            this_argset = set()
            for numruns in range(1, op_argdict["runs"]["totalmult"]+1):
                thisargs = []
                tups=gen_sumtuples(numruns, op_argdict["runs"]["totalmult"])
                for tup in tups:
                    slots=gen_kpattern_slotsets((2 * op_argdict["runs"]["loop"])+op_argdict["runs"]["totalmult"], tup, numruns)
                thisargs.append(gen_lettersets(numruns))
                #for elem in itertools.product(*thisargs): print(elem)
                this_argset |= set(itertools.product(*thisargs))
            #print(this_argset)
        else:
            #all sets of (numruns) runs with total length of runs in set < totalmult
            thisargs = []
            tups=gen_sumtuples(op_argdict["runs"]["numruns"], op_argdict["runs"]["totalmult"])
            thisargs.append(
                set(itertools.product(tup, gen_kpattern_slotsets((2 * op_argdict["runs"]["loop"])+op_argdict["runs"]["totalmult"], tup, op_argdict["runs"]["numruns"]))) for
                tup in tups)
            thisargs.append(gen_lettersets(op_argdict["runs"]["numruns"])) #what letter the runs are
            this_argset = set(itertools.product(*thisargs))
            #print(this_argset)
        op_args.append(this_argset)

    if len(op_args) > 1:
        print("Getting combined args!")
        op_args = set(itertools.product(*op_args))
    else:
        op_args= set((elem,) for elem in op_args[0])
    return op_args

def gen_argset_size(op_argdict):
    argsize=1
    if "slots" in op_argdict:
        a=2 * op_argdict["slots"]["loop"]
        k=min(op_argdict["slots"]["k"],a-1)
        n=op_argdict["slots"]["numslots"]
        argsize *= -1*((k-n+1)*(k*(n+1)-a*(n+2)-2)*binom(k+1, n))/((n+1)*(n+2))
    if "letters" in op_argdict:
        argsize*=pow(6,op_argdict["letters"]["numslots"])
    if "sumtups" in op_argdict:
        argsize *= binom(op_argdict["sumtups"]["totalmult"]-1, op_argdict["sumtups"]["numslots"]-1)
    if "rot_ind" in op_argdict:
        argsize *= 6
    return argsize
########################
# get random op_arg value
########################
def gen_random_slotcombo(keylen, k, nslots, exact=False):
    # k=1 is adjacent. if exact, max - min <= k. else, max-min == k
    if k < 1 or nslots > k: raise ValueError
    def random_slotcombo_gen():
        # generate one slot at a time
        if exact:
            this_slot = random.randrange(0, keylen - max(nslots, k))
            kbound = this_slot + k
        else:
            this_slot = random.randrange(0, keylen - nslots)
            kbound = min(this_slot + k, keylen)
        i = 1;yield this_slot
        while i < nslots - 1:
            # get the next slot.
            lastslot = this_slot
            thisbound = kbound - nslots + i + 1
            this_slot = random.randrange(lastslot + 1, thisbound + 1)
            i += 1;yield this_slot
        if exact or (this_slot + 1 == kbound):
            yield kbound
        else:
            yield random.randrange(this_slot + 1, kbound)
    return tuple(slot for slot in random_slotcombo_gen())

def gen_random_letterset(nletts):
    return tuple(alphabet[int(len(alphabet) * random.random())] for _ in range(nletts))

def gen_random_sumtuple(n_elems, target_sum):
    if n_elems > target_sum:
        print(f"cannot generate {n_elems} that sum to {target_sum}!")
    def gen_next(elems,target):
        if elems == 0: return
        elem=1+random.randrange(0,target-elems+1)
        yield elem
        yield from gen_next(elems-1,target-elem)
    return tuple(elem for elem in gen_next(n_elems, target_sum))

def get_random_argset(op_argdict):
    op_args=[]
    if "slots" in op_argdict:
        op_args.append(gen_random_slotcombo(2 * op_argdict["slots"]["loop"],
                                    op_argdict["slots"]["k"], op_argdict["slots"]["numslots"]))
    if "letters" in op_argdict:
        op_args.append(gen_random_letterset(op_argdict["letters"]["numslots"]))
    if "sumtups" in op_argdict:
        op_args.append(gen_random_sumtuple(op_argdict["sumtups"]["numslots"],op_argdict["sumtups"]["totalmult"]))
    if "rot_ind" in op_argdict:
        op_args.append(int(6*random.random()))
    return tuple(op_args)

################################
def get_mapdict(key,op_args,operation,targetsymbs={},bad_targets=None,opt='drop_bad_targets',
                argsfirst=False,no_zero_targets=False, valset=None):

    if argsfirst:
        #takes form: {src: {slot0:tgt0,slot1:tgt1,...}... etc.}
        if (opt == 'drop_bad_targets') and (bad_targets):
            fulldict = {argtup: target
                        for argtup in op_args if ((target := operation(key, *argtup)) not in bad_targets and (len(target) % 2 == 0)
                            and not (no_zero_targets and get_coeff_from_word(target,targetsymbs[int(len(target)/2)]) == 0))}
        else:
            fulldict = {argtup: target for argtup in op_args
                        if (len(target) % 2 == 0) and not (no_zero_targets and
                                                           get_coeff_from_word(target := operation(key, *argtup),
                                                                               targetsymbs[int(len(target)/2)]) == 0)}

        if (opt == 'drop_source_if_bad_targets'):
            if len(set(fulldict.values()) & set(bad_targets)) != 0: return {}
        if valset is not None:
            valset.add(target for target in fulldict.values())
    else:
        #takes form: {src: {tgt: [slot0,slot1,...slotN] etc.}
        fulldict= {}
        for argtup in op_args:
            target = operation(key, *argtup)
            if no_zero_targets:
                if (len(target) % 2 != 0): continue
                if get_coeff_from_word(target,targetsymbs[int(len(target)/2)]) == 0: continue
            if (opt == 'drop_bad_targets') and (bad_targets):
                if target in bad_targets: continue
            if (opt == 'drop_source_if_bad_targets') and (bad_targets):
                if target in bad_targets:
                    return {}

            if target in fulldict:
                fulldict[target].add(argtup)
            else:
                fulldict[target] = fastRandomSampler({argtup},inplace=True)

            if valset is not None:
                valset.add(target)
    return fastRandomSampler(fulldict,inplace=True)

def opsymb_generator(sourcesymb, targetsymbs, target_badsymb, operator, op_args, opt='drop_bad_targets', no_zero_targets=False):
    #assume we've already pruned the source symb
    valset=set()
    outdict={key:get_mapdict(key,op_args,operator,targetsymbs,target_badsymb,no_zero_targets=no_zero_targets, opt=opt, valset=valset) for key in sourcesymb}
    #print(outdict)
    print(f"Preprocessed {len(outdict)} keys in input symbol. Got {len(set(valset))} unique target keys")
    return fastRandomSampler(outdict)

def prune_opsymb(opsymb, bad_source_symb, bad_tgt_symb, drop_source_if_bad_targets=False):
    print(f"pruning: starting with {len(opsymb.keys())} source keys and {len(opsymb.values())} target keys")
    for k,v in list(opsymb.items()):
        if k in bad_source_symb: opsymb.popitem(k)
        else:
            for tgt in (v.keys() & bad_tgt_symb.keys()):
                #if tgt in opsymb[k]:
                if drop_source_if_bad_targets: opsymb.popitem(k)
                else:
                    opsymb[k].popitem(tgt)
                    if len(opsymb[k]) == 0: opsymb.popitem(k)

    print(f"After pruning, there are {len(opsymb.keys())} source keys and {len(set(opsymb.values()))} target keys")
    return opsymb

def check_key_and_get_slots(symb, loop, rel, rel_slot, format):
    #for all keys in the symb, check whether they contain the desired substring in a valid slot. If so, store the slot.
    nletter=len(list(rel.keys())[0])
    my_slot= None
    if format == "full":
        if rel_slot is not None:
            if rel_slot == -1:
                my_slot=2*loop - nletter
            else:
                my_slot = rel_slot
    elif format == "quad":
        if rel_slot is not None:
            if rel_slot == -1:
                raise ValueError
            else:
                my_slot = rel_slot + 1
            if (my_slot + nletter) > (2 * loop - 4): print("Error! bad slot!"); raise ValueError
    for symbkey in symb:
        slots = None
        if (rel_slot is not None):
            if any(check_slot(symbkey, substr, my_slot) for substr in rel):
                slots = {my_slot}
        elif format == "full":
            slots = {slot for key in rel for slot in find_all(symbkey,key)}
        elif format == "quad":
            slots = {slot+1 for key in rel for slot in find_all(symbkey[1:],key)}
        elif format == "sewmat":
            slots = {slot for key in rel for slot in find_all(symbkey[1:-1],key)}
        if not slots: continue
        yield symbkey, slots

def relsymb_generator(relnames, rels, overlaps, rel_slots, trimsymb, loop, format):
    for name, rel, rel_slot, overlap in zip(relnames, rels, rel_slots, overlaps):
        print(f'generating relsymb {name}:{datetime.datetime.now()}')
        if overlap == 0:
            yield {}
        else:
            if rel is None:
                yield trimsymb
            else:
                yield {symbkey: slots for symbkey, slots in
                       check_key_and_get_slots(trimsymb, loop, rel, rel_slot, format)}

def prune_relsymbs(relsymbs, badsymb=None):
    if badsymb is None:
        badsymb = {}
    newsymbs = [{k: relsymb[k] for k in {*relsymb} - {*badsymb}} for relsymb in relsymbs]
    return newsymbs


def tag_opinstance(instance, op_meta_args, op_tags, argslist, tagmode, is_pseudodata):
    #get the names of the operator args
    argtypes=list(op_meta_args.keys())
    #print(argtypes)
    if "slots" in op_meta_args: slotslist=argslist[argtypes.index("slots")]
    if "letters" in op_meta_args: lettslist=argslist[argtypes.index("letters")]
    if "sumtups" in op_meta_args: sumtupslist=argslist[argtypes.index("sumtups")]
    if "rot_inds" in op_meta_args: rot_indlist = argslist[argtypes.index("rot_inds")]
    my_tags = ""

    # if we're counting 'first occurrence of a', do we mean in source or target
    #if 'UNSTRIKE' in op_tags:
    #    tgt_ref = True
    #    if not tagmode: tagmode = "letters_and_slots_left"
    #else:
    #    tgt_ref = False
    #    if not tagmode: tagmode = "slots"

    if "slots" in op_meta_args and "letters" in op_meta_args:
        if not len(slotslist) == len(lettslist): return
        my_tags = op_tags
        for slot, lett in zip(slotslist, lettslist):
            my_tags = my_tags + [f'SLOT_{slot}', f'LETTER_{lett}']

    if "slots" in op_meta_args and len(op_meta_args) == 1:
        tgt_ref=False
        if not tagmode: tagmode = "slots"
        #if we are doing slots only, can do it many ways
        if tagmode == 'slots':
             my_tags=op_tags + [f'SLOT_{slot}' for slot in slotslist]
        elif (tgt_ref and len(instance["target"]) == 1) or (not tgt_ref and len(instance["source"]) == 1):
            if tgt_ref: my_key=next(iter(instance["target"]))
            else: my_key=next(iter(instance["source"]))

            if tagmode == 'letter_appearances_left':
                # i.e. STRIKE_a APP_1 means strike the first 'a' from the left
                my_tags= op_tags + [f for slot in slotslist for f in (
                                    f'LETTER_{my_key[slot]}', f'APP_{count_appearances(my_key[slot], slot)}')]
            elif tagmode == 'letter_appearances_right':
                # i.e. STRIKE_a APP_1 means strike the first 'a' from the left
                my_tags= op_tags + [f for slot in slotslist
                                                        for f in (f'LETTER_{my_key[slot]}',
                                                                  f'RAPP_{count_appearances(my_key[::-1][slot], slot)}')]
            elif tagmode == 'letters_and_slots_left':
                my_tags= op_tags + [f for slot in slotslist
                                                        for f in (f'LETTER_{my_key[slot]}', f'SLOT_{slot}')]
            elif tagmode == 'letters_and_slots_right':
                my_tags= op_tags + [f for slot in slotslist
                                                        for f in
                                                        (f'LETTER_{my_key[slot]}',
                                                         f'RSLOT_{len(my_key) - 1 - slot}')]
            elif tagmode == 'slots_right':
                my_tags= op_tags + [f'RSLOT_{len(my_key) - 1 - slot}' for slot in slotslist]
            elif tagmode == 'letters_only':
                my_tags= op_tags + [f'LETTER_{my_key[slot]}' for slot in slotslist]
            else:
                print("invalid tag mode!")
                raise ValueError
        else:
            print("invalid tag mode!")
            raise ValueError

    if is_pseudodata:
        rel_instance = {'instance': instance,
                    'tags': {'operator': my_tags,
                             'label':'PSEUDO'}}
    else:
        rel_instance = {'instance': instance,
                    'tags': {'operator': my_tags}}
    return rel_instance

def tag_rel_instance(instance,rel_tags,my_slot, is_pseudodata):
    mytags=[tag for tag in rel_tags]
    if my_slot is not None: mytags += [f'SLOT_{my_slot}']
    if is_pseudodata: mytags += ['PSEUDO']
    rel_instance = {'instance':{'source':instance},
                            'tags':{'label':mytags}}

    return rel_instance
