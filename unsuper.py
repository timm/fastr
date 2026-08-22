"""
unsuper: unsupervised recursive bi-clustering (no goals used).
Fastmap poles halve the data; contrast names the split; recurse.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license
"""
import sys; sys.dont_write_bytecode = True
from contrast import *

def distx(tbl, row1, row2): # distance between none-goal values
  d, n, p = 0, 1/BIG, the.p
  for c in tbl.cols.x:
    n, d = n+1, d + _distx(tbl.cols.all[c], row1[c], row2[c])**p
  return (d/n) ** (1/p)

def _distx(col, a, b): # helper function for distx
  if a == b == "?": return 1
  if type(col) is dict: return a != b
  a, b = norm(col, a), norm(col, b)
  if a == "?": a = 1 if b < 0.5 else 0
  if b == "?": b = 1 if a < 0.5 else 0
  return abs(a - b)

def proj(x, lohi, row): # project row onto a line lo -> hi
  a, b = x(row, lohi.lo), x(row, lohi.hi)
  return (a*a + lohi.c*lohi.c - b*b) / (2*lohi.c + 1/BIG)

def fastmap(rows, x): # sort rows on a line between 2 poles
  some = sample(rows, min(the.few, len(rows)))
  w    = choice(some)
  a    = max(some, key=lambda r: x(w, r))
  b    = max(some, key=lambda r: x(a, r))
  lohi = o(lo=a, hi=b, c=x(a, b))
  return lohi, sorted(rows, key=lambda r: proj(x, lohi, r))

def cluster(tbl, rows=None): # bi-cluster; contrast names splits
  x = lambda r1, r2: distx(tbl, r1, r2)
  def go(rows):
    t = o(rows=rows, n=len(rows))
    if len(rows) > the.stop:
      _, tmp = fastmap(rows, x)
      half = len(tmp) // 2
      if cut := contrasts(tbl, tmp[:half], tmp[half:]):
        yes = [r for r in rows if selects(cut, r)]
        no  = [r for r in rows if not selects(cut, r)]
        if 0 < len(yes) < len(rows):
          t.cut, t.yes, t.no = cut, go(yes), go(no)
    return t
  return go(tbl.rows if rows is None else rows)

def show(t, pre=None, txt="", more=lambda t: ""): # print tree
  print(f"{pre or ''}{txt}n={t.n}{more(t)}")
  if hasattr(t, "cut"):
    sub = "" if pre is None else pre + "|  "
    show(t.yes, sub, t.cut.txt + "; ", more)
    show(t.no,  sub, "!" + t.cut.txt + "; ", more)

def leaves(t): # all terminal nodes of a tree
  if hasattr(t, "cut"): return leaves(t.yes) + leaves(t.no)
  return [t]


# ---------------------------------------------------------------
def test_distx():
  "distx is symmetric, zero on self, in 0..1"
  t = Tbl(csv(the.file))
  r1, r2 = t.rows[0], t.rows[-1]
  assert distx(t, r1, r1) == 0
  assert distx(t, r1, r2) == distx(t, r2, r1)
  assert all(0 <= distx(t, r1, r) <= 1 for r in t.rows)
  print("d(r1,r1) %.2f | d(r1,r2) %.2f" %
        (distx(t, r1, r1), distx(t, r1, r2)))
  print("r1", r1, "\nr2", r2)

def test_fastmap():
  "fastmap sorts rows along a line between poles"
  t = Tbl(csv(the.file))
  x = lambda r1, r2: distx(t, r1, r2)
  lohi, rows = fastmap(t.rows[:], x)
  assert lohi.c > 0
  ps = [proj(x, lohi, r) for r in rows]
  assert ps == sorted(ps)
  print("pole1", lohi.lo, "\npole2", lohi.hi)
  print("separation %.2f; projections %.2f .. %.2f" %
        (lohi.c, ps[0], ps[-1]))

def test_cluster():
  "print bi-cluster tree; leaf counts sum to n"
  t = Tbl(csv(the.file))
  tree = cluster(t)
  assert tree.n == len(t.rows)
  assert tree.n == sum(l.n for l in leaves(tree))
  show(tree)
  print(len(leaves(tree)), "leaves")

eg = {
      "-distx": test_distx,
      "-fastmap": test_fastmap,
      "-cluster": test_cluster}

if __name__ == "__main__": main(meta(eg))
