"""
tbl: Num, Sym, Cols, Tbl: incremental column summaries.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license
"""
import sys; sys.dont_write_bytecode = True
from start import *

Sym = dict              # where to track columns of symbols
Num = lambda: (0, 0, 0) # n,mu,m2. tracks number columns

def Col(s): # Nums start with an upper case letter
  return Num() if s[0].isupper() else Sym()

def Cols(names): # names --> list[Sym|Num]; grouped into 'x,y'.
  x, y, all, klass = {}, {}, [], None
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
  if v != "?":
    if   type(i) is tuple: i = welford(*i, v, w)
    elif type(i) is dict : i[v] = i.get(v, 0) + w
    elif i.it is Cols    : i.all = [add(c, u, w)
                                    for c, u in zip(i.all, v)]
    elif i.it is XY      :
      x, y = v
      if x != "?" and y != "?":
        b = x if i.cuts is None else _bin(i.cuts, x)
        i.seen[b] = add(i.seen.get(b) or i.y(), y, w)
    elif i.it is Tbl     : i.rows += [v]; add(i.cols, v, w)
  return i

def XY(cuts=None, y=Sym): # bin x; summarize y col per bin
  "add (x,y) pairs; seen[bin] grows a Sym (or Num) of y."
  return o(it=XY, cuts=cuts, y=y, seen={})

def ranges(tbl): # per numeric x col: the.bins cuts, +/-3sd
  if not hasattr(tbl, "cuts"):
    tbl.cuts = {}
    for at in tbl.cols.x:
      col = tbl.cols.all[at]
      if type(col) is tuple:
        lo, hi = col[1] - 3*sd(col), col[1] + 3*sd(col)
        tbl.cuts[at] = [lo + k*(hi - lo)/the.bins
                        for k in range(the.bins + 1)]
  return tbl.cuts

def _bin(cuts, v): # index of the cut interval holding v
  lo, hi = cuts[0], cuts[-1]
  return max(0, min(len(cuts) - 2,
             int((v - lo)/(hi - lo + 1/BIG)*(len(cuts)-1))))

def clone(tbl, rows=()): # new empty Tbl, same column structure
  return adds(rows,
              o(it=Tbl, rows=[], cols=Cols(tbl.cols.names)))

def adds(src, i=None): # add in all values in src.
  i = Num() if i is None else i
  for v in src: i = add(i, v)
  return i

def welford(n, mu, m2, v, w): # update a Num
  n  += w
  d   = v - mu
  mu += w * d / n
  m2 += w * d * (v - mu)
  return Num() if n < 1 else (n, mu, m2)

def mids(tbl, cols=None): # return middle of columns
  return [mid(col) for col in cols or tbl.cols.all]

def mid(col): # return the middle of distributions
  return col[1] if type(col) is tuple else max(col, key=col.get)

def sd(num): # a diversity measure for Numeric columns
  n, mu, m2 = num; return 0 if n <= 1 else (m2/(n - 1))**0.5

def ent(d): # a diversity measure for Symbolic columns
  N = sum(d.values())
  return -sum(n/N*log(n/N, 2) for n in d.values() if n > 0)

def norm(col, i): # maps i ==> 0..1 via logistic cdf approx
  if i == "?": return i
  z = (i - col[1]) / (sd(col) + 1/BIG)
  return 1 / (1 + exp(-1.7 * max(-3, min(3, z))))


# ---------------------------------------------------------------
def test_num():
  "Num tracks n, mean, sd; add then subtract restores state"
  i = Num()
  for _ in range(100): i = add(i, rand())
  n, mu, m2 = i
  assert n == 100 and abs(mu - 0.5) < 0.1
  assert abs(sd(i) - 12**-0.5) < 0.05
  i = add(add(i, 0.5), 0.5, w=-1)     # add then subtract
  assert abs(i[1] - mu) < 1e-9 and i[0] == n
  print("100 rands: n %s mu %.3f sd %.3f" % (n, mu, sd(i)))

def test_sym():
  "Sym counts symbols; mid is the mode"
  i = adds("aabbbc", Sym())
  assert i["b"] == 3 and mid(i) == "b"
  print("'aabbbc' ->", i, "mid:", mid(i))

def test_tbl():
  "load the.file; rows and x,y columns found"
  t = Tbl(csv(the.file))
  assert len(t.rows) > 100 and t.cols.x and t.cols.y
  print(len(t.rows), "rows; x:", list(t.cols.x), "y:", t.cols.y)

def test_mid():
  "show middle of every column in the.file"
  t = Tbl(csv(the.file))
  m = mids(t)
  assert len(m) == len(t.cols.names)
  print(m)

def test_norm():
  "norm maps numbers to 0..1, monotonically"
  num = adds(sample(range(1000), 100), Num())
  vs = [norm(num, v) for v in range(0, 1000, 99)]
  assert all(0 <= v <= 1 for v in vs) and vs == sorted(vs)
  print(" ".join("%d->%.2f" % (x, norm(num, x))
                 for x in [0, 250, 500, 750, 999]))

def test_cdf():
  "norm(hi)-norm(lo) ~ fraction of data inside lo..hi"
  xs = [sum(rand() for _ in range(12)) - 6
        for _ in range(1000)]
  num = adds(xs, Num())
  for lo, hi in [(-1, 1), (0, 1), (-2, 2)]:
    got  = norm(num, hi) - norm(num, lo)
    want = sum(lo <= x <= hi for x in xs) / len(xs)
    print("mass %2s..%s: cdf %.2f truth %.2f" %
          (lo, hi, got, want))
    assert abs(got - want) < 0.05

def test_xy():
  "per x-range: y counts (Sym) or y mu,sd (Num)"
  t = Tbl(csv(the.file))
  at = 0                                  # x = Clndrs
  s  = t.cols.names.index("origin")       # sym y
  n  = t.cols.names.index("Lbs-")         # num y
  xy1, xy2 = XY(ranges(t)[at]), XY(ranges(t)[at], Num)
  for r in t.rows:
    add(xy1, (r[at], r[s]))
    add(xy2, (r[at], r[n]))
  assert sum(sum(d.values()) for d in xy1.seen.values()) \
         == len(t.rows)
  for b in sorted(xy1.seen):
    print("bin %2s origin %-22s Lbs mu %5.0f sd %4.0f" %
          (b, xy1.seen[b], xy2.seen[b][1], sd(xy2.seen[b])))

eg = {
      "-num": test_num,
      "-sym": test_sym,
      "-tbl": test_tbl,
      "-mid": test_mid,
      "-norm": test_norm,
      "-cdf": test_cdf,
      "-xy": test_xy}

if __name__ == "__main__": main(meta(eg))
