import os,re
import numpy as np
import torch
import sympy
from fractions import Fraction

#some classes to hold the data in different formats.
#symb is an overload of dict with some elementwise operators on values,
#sumlist is an overload of list with elementwise sum and multiplication operations

class Symb(dict):
    def valmult(self,m,d1):
        return {k:(m*v) for k,v in d1.items()}
    
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
    
    def __getitem__(self,key):
        if key in self: return super().__getitem__(key)
        else: return 0
    
class sumlist():
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

    def __rsub__(self, other):
        if isinstance(otherlist,sumlist): os=otherlist.list
        else: os=otherlist
        
        outlist=[a_i - b_i for a_i, b_i in zip(otherlist, os)]
        return sumlist(outlist)
    
    def __mul__(self,const):
        return [const*elem for elem in self.list]

    def __rmul__(self,const):
        return [const*elem for elem in self.list]
    
    def __getitem__(self,key):
        if key in self: return super().__getitem__(key)
        else: return 0
