"""
Core data structures for symbol manipulation.

Symb: dict subclass with arithmetic operator overloads (+, -, *, /, &, |)
for elementwise operations on symbol coefficients. Supports in-place
add_small/sub_small for incremental updates without creating new dicts.

sumlist: list wrapper with elementwise addition/subtraction/multiplication,
used for lists of Symb objects (e.g. coproduct operations in common_dev).

fastRandomSampler: O(1) random sampling + O(1) lookup data structure.
Wraps a dict or set with a bidirectional key<->int mapping so that random
key selection and key removal are both constant-time. Critical for efficient
training data generation from symbols with millions of terms.
"""

import copy
import random
import argparse
from fractions import Fraction


class Symb(dict):
    """Dict subclass representing a symbol (linear combination of letter-string keys).

    Keys are letter strings (e.g. 'aabbcd'), values are integer coefficients.
    Arithmetic operators act on coefficients: (s1 + s2)[k] = s1[k] + s2[k].
    Zero-valued entries are automatically dropped from addition/subtraction results.
    """
    def dict(self):
        return {k:v for k,v in self.items()}

    def int_if(self,x):
        return int(x) if float(x).is_integer() else x

    def valmult(self,m,d1):
        return {k:self.int_if(m*v) for k,v in d1.items()}

    def valdiv(self,m,d1):
        return {k:self.int_if(v/m) for k,v in d1.items()}

    def intcast(self):
        return Symb(self.dict())

    def dictmerge(self,d1,d2):
        def getval(k):
            if k in d1 and k in d2: val= d1[k]+d2[k]
            elif k in d1: val= d1[k]
            elif k in d2: val= d2[k]
            else: val = None
            return val

        return{k:getval(k) for k in d1|d2 if getval(k) != 0}
    
    def dictdiff(self,d1,d2):
        def getval(k):
            if k in d1 and k in d2: val= d1[k]-d2[k]
            elif k in d1: val= d1[k]
            elif k in d2: val= -d2[k]
            else: val = None
            return val
    
        return{k:getval(k) for k in d1|d2 if getval(k) != 0}

    def add_small(self, othersymb):
        if isinstance(othersymb,Symb): os=othersymb
        else: os=Symb(othersymb)

        def getval(k):
            if k in self and k in othersymb:
                val = self[k] + othersymb[k]
            elif k in self:
                val = self[k]
            elif k in othersymb:
                val = othersymb[k]
            else:
                val = None
            return val

        for k in othersymb: self[k] = getval(k)
        return self

    def sub_small(self, othersymb):
        #subtract the othersymb from self
        if isinstance(othersymb, Symb):
            os = othersymb
        else:
            os = Symb(othersymb)

        def getval(k):
            if k in self and k in othersymb:
                val = self[k] - othersymb[k]
            elif k in self:
                val = self[k]
            elif k in othersymb:
                val = -othersymb[k]
            else:
                val = None
            return val

        for k in othersymb: self[k] = getval(k)
        return self

    def __add__(self, othersymb):
        if isinstance(othersymb,Symb): os=othersymb
        else: os=Symb(othersymb)
        s=self.dictmerge(self,os)
        #print(Symb(s))
        #print(s,os.mydict,self.mydict)
        return Symb(s)
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, othersymb):
        return Symb(self.dictdiff(self,othersymb))

    def __rsub__(self, othersymb):
        return Symb(self.dictdiff(othersymb,self))

    def __mul__(self,const):
        return Symb(self.valmult(const,self))
    
    def __rmul__(self,const):
        return Symb(self.valmult(const,self))

    def __truediv__(self,const):
        return Symb(self.valdiv(const,self))

    def __and__(self, othersymb):
        if isinstance(othersymb,Symb):
            return Symb(self.dict() & othersymb.dict())
        else:
            return Symb(self.dict() & othersymb)

    def __rand__(self, othersymb):
        return self.__and__(othersymb)

    def __or__(self, othersymb):
        if isinstance(othersymb,Symb):
            return Symb(self.dict() | othersymb.dict())
        else:
            return Symb(self.dict() | othersymb)

    def __ror__(self, othersymb):
        return self.__or__(othersymb)


class sumlist():
    """List wrapper with elementwise Symb arithmetic.

    Used for lists of Symb dicts where (list1 + list2)[i] = list1[i] + list2[i].
    Used by coproduct operations as collections of per-basis-element constraints.
    """
    def __init__(self,mylist):
        self.list=mylist
        
    def __add__(self, otherlist):
        if isinstance(otherlist,sumlist): os=otherlist.list
        else: os=otherlist
        
        outlist=[a_i + b_i for a_i, b_i in zip(self.list, os)]
        return sumlist(outlist)

    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, otherlist):
        if isinstance(otherlist,sumlist): os=otherlist.list
        else: os=otherlist
        
        outlist=[a_i - b_i for a_i, b_i in zip(self.list, os)]
        return sumlist(outlist)

    def __rsub__(self, otherlist):
        if isinstance(otherlist,sumlist): os=otherlist.list
        else: os=otherlist

        outlist=[a_i - b_i for a_i, b_i in zip(os, self.list)]
        return sumlist(outlist)
    
    def __mul__(self,const):
        return [const*elem for elem in self.list]

    def __rmul__(self,const):
        return [const*elem for elem in self.list]
    
    def __getitem__(self,key):
        return self.list[key]

class fastRandomSampler(object):
    """O(1) random sampling + O(1) lookup data structure for dicts and sets.

    Maintains a parallel list (keylist) and reverse map (key_to_int) alongside
    the underlying dict/set. To pop a random element:
      1. Draw a random int i, get key = keylist[i]
      2. Swap keylist[i] with the last element, pop the last
      3. Remove key from the dict/set (O(1) lookup)

    Optional countdict tracks multiplicity -- keys can be sampled multiple times
    before being removed (used for scramble/overlap tracking).

    Args:
        init_elem: dict or set to wrap.
        countdict: Optional multiplicity counts per key.
        inplace: If True, mutates init_elem directly (faster but destructive).
    """

    def __init__(self, init_elem, countdict={}, inplace=False):
        # this sampling struct works for dicts and sets, so just flag which it is
        if (not isinstance(init_elem, dict) and not isinstance(init_elem, set)): raise TypeError
        self.is_dict = isinstance(init_elem, dict)

        #the object to sample from. if inplace, edits the original dict,
        # otherwise, edits a copy. inplace is faster, but more dangerous.
        if inplace: self.mystruct = init_elem
        else: self.mystruct = init_elem.copy()

        # optional counter, if items have multiplicity. Used for scramble
        self.countdict = countdict

        # Create the key-to-int maps from the dictionary
        if len(self.mystruct) == 0:
            self.keylist, self.key_to_int = [],{}
        else:

            #self.keylist =
            self.keylist, self.key_to_int = [([*tup] if i == 0 else dict(tup))
                                                           for i,tup in enumerate(zip(*((k, (k, v))
                                                           for v, k in enumerate(self.mystruct))))]
    def __getitem__(self, item):
        if item in self.mystruct:
            return (self.mystruct[item] if self.is_dict else item)
        else:
            return None

    def __contains__(self, item):
        return item in self.mystruct

    def __len__(self):
        return len(self.mystruct)

    def __repr__(self):
        return str(self.mystruct)

    def __str__(self):
        return str(self.mystruct)

    def copy(self):
        return copy.copy(self)

    def keys(self):
        return self.mystruct.keys() if self.is_dict else self.mystruct

    def items(self):
        return self.mystruct.items() if self.is_dict else self.mystruct

    def values(self):
        return self.mystruct.values() if self.is_dict else None

    def add(self, key, value=None):  # O(1)
        # Add key-value pair (no extra work needed for simply changing the value)
        new_int = len(self.mystruct)
        if self.is_dict:
            self.mystruct[key] = value
        else:
            self.mystruct.add(key)

        self.key_to_int[key] = new_int
        self.keylist.append(key)

    def popitem(self, key):  # O(1)

        position = self.key_to_int.pop(key)
        last_item = self.keylist.pop()

        if position != len(self.keylist):
            self.keylist[position] = last_item
            self.key_to_int[last_item] = position
        if self.is_dict:
            return key, self.mystruct.pop(key)
        else:
            self.mystruct.remove(key)
            return key

    def remove(self, key):# O(1)
        if key not in self.key_to_int: return
        self.popitem(key)
        return

    def random_key(self):  # O(1)
    # Select a random key from the dictionary using the int_to_key map
        return self.keylist[int(len(self.mystruct) * random.random())]

    def remove_random(self):  # O(1)
        # Randomly remove a key from the dictionary via the bidirectional maps
        key = self.random_key()
        self.remove(key)

    def pop_random(self):  # O(1)
        # Randomly pop a key from the dictionary via the bidirectional maps
        try:
            key = self.random_key()
            if not self.countdict:
                # if we're not counting, just pop it
                return self.popitem(key)
            else:
                # countdict lets us pop keys multiple times (if keys have multiplicity)
                if self.countdict[key] == 1:
                    # If we're on the last time, pop it
                    self.countdict.pop(key, None)
                    return self.popitem(key)
                else:
                    # otherwise, decrement its counter- we've seen it
                    self.countdict[key] -= 1
                    if self.is_dict:
                        return key, self.mystruct[key]
                    else:
                        return key
        except IndexError:
            print("Error, symbol is exhausted!")

    def pop_random_gen(self, num_to_pop):  # O(1)
        for i in range(num_to_pop):
            yield self.pop_random()

    def pop_inst_gen(self, subdict_size, num_to_gen):
        # Randomly pop instances from the symb
        for i in range(num_to_gen):
            if self.is_dict:
                yield {k: v for k, v in self.pop_random_gen(subdict_size)}
            else:
                yield {k for k in self.pop_random_gen(subdict_size)}

FALSY_STRINGS = {'off', 'false', '0'}
TRUTHY_STRINGS = {'on', 'true', '1'}
def bool_flag(s):
    """
    Parse boolean arguments from the command line.
    """
    if s.lower() in FALSY_STRINGS:
        return False
    elif s.lower() in TRUTHY_STRINGS:
        return True
    else:
        raise argparse.ArgumentTypeError("Invalid value for a boolean flag!")

