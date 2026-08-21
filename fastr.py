"""
shortr: tiny multi-obkective explainable AI.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license

Options:
 
  -h             show help
  --p=2          minkowski coeffeience
  --stop=4       stopping rule for recursive tree generation
  --few=256      sub-sample size for pole finding
  --k=5          nearest neighbors used in a leaf
  --check=5      optimization: how many top picks to evaluate
  --seed=1234    random number genration
  --file=/Users/timm/gits/moot/optimize/misc/auto93.csv
 """ 
import re, sys; sys.dont_write_bytecode = True
from math import exp,cos,log
from types import SimpleNamespace as o
from random import choice,sample,seed,shuffle,random as rand

BIG=1e32

# ---------------------------------------------------------------
def atom(s): # string --> atom
  try: return int(s)
  except ValueError:
    try: return float(s)
    except ValueError: return s.strip()

def csv(f): # csv file strings --> rows of atoms
  with open(f, encoding="utf-8") as fs: 
    for s in fs:
      if s := s.strip(): yield [atom(x) for x in s.split(",")]

Sym=dict            # where to track columns of symbols
Num=lambda: (0,0,0) # n mu m2. where to tracjkcolumns of numbers

def Col(s): # Nums start with an upper case letter
  return Num() if s[0].isupper() else Sym()

def Cols(names): # names --> list[Sym|Num]; grouped into 'x,y'.
  x,y,all,klass = {},{},[],None
  for at, name in enumerate(names):
    all += [Col(name)]
    if name[-1] == "X": continue
    elif name[-1] == "!": klass = at
    else: (y if name[-1] in "+-" else x)[at] = name[-1] == "+"
  return o(it=Cols, all=all, x=x, y=y, klass=klass, names=names)

def Tbl(src): # where to track rows and columns
  src = iter(src)
  return adds(src, o(it=Tbl, rows=[], cols=Cols(next(src))))

def add(i, v, w=1): # add one value
  if v!="?":
    if   type(i) is tuple: i = welford(*i, v, w)
    elif type(i) is dict : i[v] = i.get(v,0) + w
    elif i.it is Cols    : i.all = [add(c,v1,w) for c,v1 in zip(i.all,v)]
    elif i.it is Tbl     : i.rows += [v]; add(i.cols,v, w)
  return i

def clone(tbl, rows=()): # new empty Tbl, same column structure
  return adds(rows, o(it=Tbl, rows=[], cols=Cols(tbl.cols.names)))

def adds(src, i=None): # add in all values in src.
  i = Num() if i is None else i
  for v in src: i = add(i, v)
  return i

def welford(n,mu,m2, v, w):  # update a Num
  n += w
  if n < 1: return Num()
  d   = v - mu
  mu += w * d / n
  m2 += w * d * (v - mu)
  return (n, mu, m2)

def mids(tbl, cols=None): # return middle of columns
  return [mid(col) for col in cols or tbl.cols.all]

def mid(col):  # return the middle of distributions
  return col[1] if type(col) is tuple else max(col, key=col.get) 

def sd(num): # a diversity measure for Numeric columns
  n,mu,m2 = num; return 0 if n <= 1 else (m2/(n - 1))**0.5 

def ent(d): # a diversity measure for Symbolic columns
  N = sum(d.values())
  return -sum(n/N*log(n/N,2) for n in d.values() if n > 0)

# ---------------------------------------------------------------
def disty(tbl, row): # distance of goals to the reference optimum
  d,n,p = 0, 1/BIG, the.p
  for c,w in tbl.cols.y.items():
    v = norm(tbl.cols.all[c], row[c])
    if v != "?": n,d = n + 1, d + abs(v - w)**p
  return (d/n) ** (1/p)

def distx(tbl, row1,row2): # distance between none-goal balues
  d,n,p = 0, 1/BIG, the.p
  for c in tbl.cols.x:
    n,d = n + 1, d + distxq(tbl.cols.all[c], row1[c],row2[c])**p
  return (d/n) ** (1/p)

def distxq(col,a,b):
  if a==b=="?": return 1
  if type(col) is dict: return a != b
  a,b = norm(col,a), norm(col,b)
  if a=="?": a = 1 if b < 0.5 else 0
  if b=="?": b = 1 if a < 0.5 else 0
  return abs(a-b)

def norm(col, i): # maps i ==> 0..1 via logistic cdf approx
  if i=="?": return i
  z = (i - col[1]) / (sd(col) + 1/BIG)
  return 1 / (1 + exp(-1.7 * max(-3, min(3, z))))

# ---------------------------------------------------------------
def confuse(pairs): # [(got,want)] --> per-class acc,pd,pf,prec
  out = {k: o(label=k, tp=0, fp=0, fn=0, tn=0)
         for g, w in pairs for k in (g, w)}
  for got, want in pairs:
    for k, c in out.items():
      if   k == want: c.tp += got == want; c.fn += got != want
      elif k == got : c.fp += 1
      else          : c.tn += 1
  for c in out.values():
    c.acc  = (c.tp + c.tn) / (c.tp+c.fn+c.fp+c.tn + 1/BIG)
    c.pd   = c.tp / (c.tp + c.fn + 1/BIG)  # recall
    c.pf   = c.fp / (c.fp + c.tn + 1/BIG)  # false alarm
    c.prec = c.tp / (c.tp + c.fp + 1/BIG)
  return out

def xval(rows, n=5, m=5): # m repeats of n-fold cross-val
  rows = rows[:]
  for _ in range(m):
    shuffle(rows)
    for j in range(n):
      yield ([r for i, r in enumerate(rows) if i % n != j],
             rows[j::n])                                  

# ---------------------------------------------------------------
def proj(x, lohi, row): # project x onto a line a -> b
  a,b = x(row, lohi.lo), x(row, lohi.hi)
  return (a*a + lohi.c*lohi.c - b*b) / (2*lohi.c + 1/BIG)

def fastmap(rows, x, y):
  some = sample(rows, min(the.few, len(rows)))
  w    = choice(some)
  a    = max(some, key=lambda r: x(w, r))
  b    = max(some, key=lambda r: x(a, r))
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

def predict(tbl, tree, row, k=None): # klass of k nearest in leaf
  rows = sorted(leaf(tbl, tree, row).rows,
                key=lambda r: distx(tbl, row, r))[:k or the.k]
  at = tbl.cols.klass
  return mid(adds((r[at] for r in rows), Col(tbl.cols.names[at])))

def leaf(tbl, t, row):
  x = lambda r1, r2: distx(tbl, r1, r2)
  while hasattr(t, "mid"):
    ok = proj(x, t.lohi, row) <= t.mid
    t = t.goods if ok else t.bads
  return t
 
# ---------------------------------------------------------------
def test_h()     : print(__doc__)
def test_the()   : print(the)
def test__seed(s): the.seed = s
def test__file(f): the.file = f

def test_atom():
  assert atom("2") == 2 and atom("2.1") == 2.1
  assert atom(" a ") == "a"

def test_num():
  i = Num()
  for _ in range(100): i = add(i, rand())
  n,mu,m2 = i
  assert n == 100 and abs(mu - 0.5) < 0.1
  assert abs(sd(i) - 12**-0.5) < 0.05
  i = add(add(i, 0.5), 0.5, w=-1)     # add then subtract
  assert abs(i[1] - mu) < 1e-9 and i[0] == n

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
  print(mids(Tbl(csv(the.file))))

def test_norm():
  num = adds(sample(range(1000), 100), Num())
  vs = [norm(num, v) for v in range(0, 1000, 99)]
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

def test_confuse():
  pairs = [("a","a")]*4 + [("b","a")] + \
          [("b","b")]*2 + [("a","b")]
  a = confuse(pairs)["a"]
  assert a.tp == 4 and a.fn == 1 and a.fp == 1 and a.tn == 2
  assert a.pd == 0.8 and abs(a.prec - 0.8) < 1e-9
  for c in confuse(pairs).values():
    print(c.label, "acc %.2f pd %.2f pf %.2f prec %.2f" %
          (c.acc, c.pd, c.pf, c.prec))

def test_xval():
  for f in ["diabetes", "soybean"]:
    print("\n##", f)
    t = Tbl(csv(f"/Users/timm/gits/moot/classify/{f}.csv"))
    k, pairs = t.cols.klass, []
    for train, test in xval(t.rows, n=3, m=10):
      tr   = clone(t, train)
      tree = fastmapr(tr, tr.rows, stop=128)
      for row in test:
        pairs += [(predict(tr, tree, row), row[k])]
    for c in sorted(confuse(pairs).values(), key=lambda c: -c.pd):
      print("%25s acc %.2f pd %.2f pf %.2f prec %.2f" %
            (c.label, c.acc, c.pd, c.pf, c.prec))

def test_opt():
  t = Tbl(csv(the.file))
  y = lambda tbl, r: disty(tbl, r)
  best, rnd = Num(), Num()
  for train, test in xval(t.rows, n=3, m=10):
    tr   = clone(t, train)
    tree = fastmapr(tr, tr.rows)
    score = lambda r: mid(adds((y(tr, r1)
                    for r1 in leaf(tr, tree, r).rows), Num()))
    picks = sorted(test, key=score)[:the.check]
    best  = add(best, min(y(tr, r) for r in picks))
    rnd   = add(rnd,  min(y(tr, r) for r in sample(test, the.check)))
  lo = min(y(t, r) for r in t.rows)
  print("true best %.2f | best of %s picked %.2f | random %.2f" %
        (lo, the.check, mid(best), mid(rnd)))

def test_all():
  for s, f in list(globals().items()):
    if re.match(r"test_(?!_|all)", s): print("\n#", s); run(f)

# ---------------------------------------------------------------
def run(f, *v):
  try: seed(the.seed); f(*v)
  except Exception: import traceback; traceback.print_exc()

the=o(**{k:atom(v) for k,v in re.findall(r"(\w+)=(\S+)",__doc__)})
seed(the.seed)

def main():
  for s, v in zip(sys.argv, sys.argv[1:] + [None]):
    if f := globals().get(f"test{s.replace('-', '_')}"):
      run(f, atom(v)) if s[1] == "-" else run(f)

if __name__ == "__main__": main()

