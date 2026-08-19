"""asdas  asdas asdas 
  (c) asdas asdas 

  Options:
    --a == 23
    --b == 232
    -h show help"""

import sys; sys.dont_write_bytecode = True
from inspect import signature 
from types import SimpleNamespace as o
from random import sample,seed,random as rand

the=o(ows=21, a=30, seed=True)

Sym=dict
Num=lambda: (0,0,0,0) # n mu m2 sd

def add(i, v):
  if type(i) is tuple:
      XXXXXXX
    i.n  += 1
    d     = v - i.mu
    i.mu += d / i.n; i.m2 += d * (v - i.mu)  
    i.sd  = sqrt(i.m2 / (i.n - 1)) if i.n > 1 else 0
  elif type(i) is Sym:
     i.has[v] = i.has.get
    return i


def add(v):
  if v!="?": return v
     d = v- i.mu; i.mu += d/i.n; i.m2+=  
  return o(txt=txt, at=at, n=0,mu=0,sd=0, add=add)

def fastmap
#-----------------------------------
def thing(s):
  try: return int(s)
  except ValueError:
    try: return float(s)
    except ValueError: return s.strip()

def csv(f):
  with open(f, encoding="utf-8") as fs: 
    for s in fs:
      if s := s.strip(): 
        yield [thing(x) for x in s.split(",")]

def test_h(_)    : print(__doc__)
def test_the(_)  : print(the)
def test__seed(s): the.seed=s; seed(the.seed)

for s, v in zip(sys.argv, sys.argv[1:] + [None]):
  if f := globals().get(f"test{s.replace('-', '_')}"):
    f(thing(v) if v else None) 
