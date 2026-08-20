"""
shortr: tiny multi-obkective explainable AI.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license

Options:
 
  -h             show help
  --p=2          minkowski coeffeience
  --stop=4       stopping rule for recursive tree generation
  --seed=1234    random number genration
  --file=/Users/timm/gits/moot/optimize/misc/auto93.csv""" 

import re, sys; sys.dont_write_bytecode = True
from inspect import signature 
from types import SimpleNamespace as o
from random import choice,sample,seed,random as rand
BIG=1e32

# ---------------------------------------------------------------
Sym=dict

Num=lambda: (0,BIG,0,-BIG) # n lo mu hi

def Tbl(src):
  src = iter(src)
  return adds(src, o(rows=[], cols=Cols(next(src))))

def Cols(names): # turn a list of names into Syms or Nums
  all,x,y=[],{},{}
  for at, s in enumerate(names):
    all += [(Num if s[0].isupper() else Sym)()]
    if s[-1] == "X": continue
    (y if s[-1] in "+-" else x)[at] = s[-1] == "+"
  return o(all=all, x=x, y=y, names=names)

def adds(src, i):
  i = i or Num()
  for v in src: i = add(i, v)
  return i

def add(i, v):
  if v=="?": return i
  if type(i) is dict: i[v] = i.get(v,0) + 1
  elif type(i) is tuple:
    n,lo,mu,hi = i
    n += 1
    i = (n, min(lo,v), mu+(v-mu)/n, max(hi,v))
  else: 
    i.rows += [v]
    i.cols.all = [add(c,v1) for c,v1 in zip(i.cols.all,v)]
  return i

def mid(i):
  if type(i) is dict: return max(i, key=i.get)
  elif type(i) is tuple: return i[2]
  else: return [mid(col) for col in i.cols.all]

# def spans(best, rest, more=()):
#   cuts = sorted({best.lo, best.mu, best.hi, rest.lo, rest.mu, rest.hi, *more})
#   return [(cuts[i], cuts[j]) for i in range(len(cuts)) for j in range(i+1, len(cuts))]
#
# pick = max(spans(best, rest), key=lambda ab: score(ab, best, rest))
#

# ---------------------------------------------------------------
def disty(tbl, row):
  d,n,p = 0, 1/BIG, the.p
  for c,w in tbl.cols.y.items():
    v = norm(row[c], *tbl.cols.all[c][1:])
    if v != "?": n,d = n + 1, d + abs(v - w)**p
  return (d/n) ** (1/p)

def distx(tbl, row1,row2):
  d,n,p = 0, 1/BIG, the.p
  for c in tbl.cols.x:
    n,d = n + 1, d + distx1(tbl.cols.all[c], row1[c],row2[c])**p
  return (d/n) ** (1/p)

def distx1(col,a,b):
  if a==b=="?": return 1
  if type(col) is dict: return a != b
  a,b = norm(a,*col[1:]), norm(b,*col[1:])
  if a=="?": a = 1 if b < 0.5 else 0
  if b=="?": b = 1 if a < 0.5 else 0
  return abs(a-b)

def norm(v, lo, mu, hi): # maps v ==> 0..1
  if v=="?": return v
  if v <= lo: return 0
  if v >= hi: return 1
  w = hi - lo
  if v <= mu: return (v-lo)**2 / (w*(mu-lo))
  return 1 - (hi-v)**2 / (w*(hi-mu))

# ---------------------------------------------------------------
def proj(x, lohi, row): 
  a,b = x(row, lohi.lo), x(row, lohi.hi)
  return (a*a + lohi.c*lohi.c - b*b) / (2*lohi.c + 1/BIG)

def fastmap(rows, x, y):
  w = choice(rows)
  a = max(rows, key=lambda r: x(w, r))
  b = max(rows, key=lambda r: x(a, r))
  if y(b) < y(a): a, b = b, a
  lohi = o(lo=a, hi=b, c=x(a, b))
  return lohi, sorted(rows, key=lambda r: proj(x, lohi, r))

def fastmapr(tbl, rows, stop=None):
  x = lambda r1, r2: distx(tbl, r1, r2)
  y = lambda r: disty(tbl, r)
  stop = stop or the.stop
  def go(rows):
    t = o(rows=rows)
    if len(rows) > stop:
      lohi, t.rows = fastmap(rows, x, y)
      n = len(t.rows) // 2
      t.lohi, t.mid = lohi, proj(x, lohi, t.rows[n])
      t.goods, t.bads = go(t.rows[:n]), go(t.rows[n:])
    return t
  return go(rows)

def leaf(tbl, t, row):
  x = lambda r1, r2: distx(tbl, r1, r2)
  while hasattr(t, "mid"):
    ok = proj(x, t.lohi, row) <= t.mid
    t = t.goods if ok else t.bads
  return t
 
# ---------------------------------------------------------------
def it(s):
  try: return int(s)
  except ValueError:
    try: return float(s)
    except ValueError: return s.strip()

def csv(f):
  with open(f, encoding="utf-8") as fs: 
    for s in fs:
      if s := s.strip(): yield [it(x) for x in s.split(",")]

# ---------------------------------------------------------------
def test_h()     : print(__doc__)
def test_the()   : print(the)
def test_seed(s) : the.seed = s
def test_file(f) : the.file = f

def test_it():
  assert it("2") == 2 and it("2.1") == 2.1
  assert it(" a ") == "a"

def test_num():
  i = Num()
  for _ in range(100): i = add(i, rand())
  n,lo,mu,hi = i
  assert n == 100 and 0 <= lo <= mu <= hi <= 1
  assert abs(mu - 0.5) < 0.1

def test_sym():
  i = adds("aabbbc", Sym())
  assert i["b"] == 3 and mid(i) == "b"

def test_csv():
  assert len(list(csv(the.file))) > 100

def test_tbl():
  t = Tbl(csv(the.file))
  assert len(t.rows) > 100 and t.cols.x and t.cols.y
  print(len(t.rows), "rows; x:", list(t.cols.x), "y:", t.cols.y)

def test_mid():
  print(mid(Tbl(csv(the.file))))

def test_norm():
  i = adds(sample(range(1000), 100), Num())
  vs = [norm(v, *i[1:]) for v in range(0, 1000, 99)]
  assert all(0 <= v <= 1 for v in vs) and vs == sorted(vs)

def test_distx():
  t = Tbl(csv(the.file))
  r1, r2 = t.rows[0], t.rows[-1]
  assert distx(t, r1, r1) == 0
  assert distx(t, r1, r2) == distx(t, r2, r1)
  assert all(0 <= distx(t, r1, r) <= 1 for r in t.rows)

def test_disty():
  t = Tbl(csv(the.file))
  lst = sorted(t.rows, key=lambda r: disty(t, r))
  n=len(lst)
  print(t.cols.names)
  for j in [0,1,2,n-3,n-2,n-1]: print(t.rows[j])

def test_fastmap():
  t = Tbl(csv(the.file))
  x = lambda r1, r2: distx(t, r1, r2)
  y = lambda r: disty(t, r)
  lohi, rows = fastmap(t.rows[:], x, y)
  assert y(lohi.lo) <= y(lohi.hi) and lohi.c > 0
  ps = [proj(x, lohi, r) for r in rows]
  assert ps == sorted(ps)

def test_tree():
  t = Tbl(csv(the.file))
  def leaves(n): 
    return leaves(n.goods) + leaves(n.bads) if hasattr(n,"mid") else [n]
  tree = fastmapr(t, t.rows[:])
  ls = leaves(tree)
  assert all(len(l.rows) <= the.stop for l in ls)
  print(len(ls), "leaves")

def test_leaf():
  t = Tbl(csv(the.file))
  tree = fastmapr(t, t.rows[:])
  l = leaf(t, tree, t.rows[0])
  assert len(l.rows) <= the.stop
  print("row 0 lands with", len(l.rows), "rows")

def test_all():
  for s, f in list(globals().items()):
    if s.startswith("test_") and f is not test_all \
       and not signature(f).parameters:
      print("\n#", s); run(f)

# ---------------------------------------------------------------
def run(f, v=None):
  try: seed(the.seed); f(v) if signature(f).parameters else f()
  except Exception: import traceback; traceback.print_exc()

the=o(**{k:it(v) for k,v in re.findall(r"(\w+)=(\S+)",__doc__)})
seed(the.seed)

if __name__ == "__main__":
  for s, v in zip(sys.argv, sys.argv[1:] + [None]):
    if f := globals().get(f"test{s.replace('-', '_')}"):
      run(f, it(v) if v else None) 
