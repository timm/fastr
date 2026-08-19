"""asdas  asdas asdas 
  (c) asdas asdas 

  Options:
    --a == 23
    --b == 232
    -h show help"""

import sys; sys.dont_write_bytecode = True
from inspect import signature 
from types import SimpleNamespace as o
from random import choice,sample,seed,random as rand

the=o(ows=21, a=30, seed=True)
BIG=1e32

# --------------------------------------------------------------------
Sym=dict

Num=lambda: (BIG,0,-BIG) # lo mu hi

def Tbl(src):
  src = iter(src)
  return adds(src, o(rows=[], cols=Cols(next(src))))

def Cols(names): # turn a list of names into Syms or Nums
  all,x,y=[],{},{}
  for at, s in enumerate(names):
    cols[at] = [(Num if s[0].isupper() else Sym)()]
    if s[-1] == "X": continue
    (y if s[-1] in "+-" else x)[at] = s[-1] == "+"
  return o(all=all, x=x, y=y)

def add(i, v):
  if v=="?": return i
  if type(i) is dict: i[v] = i.get(v,0) + 1
  elif type(i) is tuple:
    lo,mu,hi = i
    i = (min(lo,v), mu+(v-mu)/n, max(hi,v))
  else: 
    i.rows += [v]
    i.cols.all = [add[c,v1) for c,v1 in zip(i.cols.all,v)]
  return i

def mid(i):
  if type(i) is dict: return max(i, key=i.get)
  elif type(i) is tuple: return i[1]
  else return [mid(col) for col in i.cols.all.values()]


# def spans(best, rest, more=()):
#   cuts = sorted({best.lo, best.mu, best.hi, rest.lo, rest.mu, rest.hi, *more})
#   return [(cuts[i], cuts[j]) for i in range(len(cuts)) for j in range(i+1, len(cuts))]
#
# pick = max(spans(best, rest), key=lambda ab: score(ab, best, rest))
#

# --------------------------------------------------------------------
def disty(tbl, row,):
  d,n,p = 0, 1/BIG, the.p
  for c,w in tble.cols.y.items():
    n,d = n + 1, d+ (norm(tbl.cols.all[c],row[c]) - w)**p
  return (d/n) ** (1/p)

def distx(tbl, row1,row2):
  d,n,p = 0, 1/BIG, the.p
  for c in tbl.cols.x:
    n,d = n + 1, d + distx1(tbl.cols.all[c], row1[c],row2[c])**p
  return (d/n) ** (1/p)

def distx1(col,a,b):
  if a==b=="?": return 1
  a,b = norm(col,a), norm(col,b)
  if a=="?": a = 1 if b < 0.5 else 0
  if b=="?": a = 1 if a < 0.5 else 0
  return abs(a-b)

def norm(v, lo, mu, hi): # maps v ==> 0..1
  if v=="?": return v
  if v <= lo: return 0
  if v >= hi: return 1
  w = hi - lo
  if v <= mu: return (v-lo)**2 / (w*(mu-lo))
  return 1 - (hi-v)**2 / (w*(hi-mu))

# --------------------------------------------------------------------

def cosine(a, b, c): return (a*a + c*c - b*b) / (2*c + 1/BIG)

def proj(tbl,a,b,c, row):
  return cosine(distx(tbl, row, a), distx(tbl, row, b), c)

def fastmapr(tbl, rows, stop=4):
  if len(rows) <= stop: return o(rows=rows)
  w = choice(rows)
  a = max(rows, key=lambda r: distx(tbl, w, r))
  b = max(rows, key=lambda r: distx(tbl, a, r))
  if disty(tbl, b) < disty(tbl, a): a, b = b, a
  rows.sort(key=lambda r: proj(tbl, a,b,c, r))
  n = len(rows) // 2
  return o(rows=rows, a=a, b=b, c=c,
           mid   = proj(tbl, here, rows[n]),
           goods = fastmapr(tbl, rows[:n], stop),
           bads  = fastmapr(tbl, rows[n:], stop))

def leaf(tbl, node, row):
  while hasattr(node, "mid"):
    ok   = proj(tbl, node, row) <= node.mid
    node = node.goods if ok else node.bads
  return node
  
# --------------------------------------------------------------------
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

#-----------------------------------
def test_h(_)    : print(__doc__)
def test_the(_)  : print(the)
def test__seed(s): the.seed=s; seed(the.seed)

for s, v in zip(sys.argv, sys.argv[1:] + [None]):
  if f := globals().get(f"test{s.replace('-', '_')}"):
    f(thing(v) if v else None) 
