"""
super: supervised recursive splitting for classify and regress.
Each split minimizes the spread of the goals in two halves.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license
"""
import sys; sys.dont_write_bytecode = True
from unsuper import *

def disty(tbl, row): # distance of goals to the reference optimum
  d, n, p = 0, 1/BIG, the.p
  for c, w in tbl.cols.y.items():
    v = norm(tbl.cols.all[c], row[c])
    if v != "?": n, d = n + 1, d + abs(v - w)**p
  return (d/n) ** (1/p)

def yinfo(tbl): # goal accessor, summary maker, spread measure
  if tbl.cols.y: return (lambda r: disty(tbl, r)), Num, sd
  at = tbl.cols.klass
  return (lambda r: r[at]), Sym, ent

def split(tbl, rows): # col+cut minimizing spread of 2 halves
  Y, New, div = yinfo(tbl)
  best = None
  for at in tbl.cols.x:
    if type(tbl.cols.all[at]) is not tuple: continue
    have = sorted((r for r in rows if r[at] != "?"),
                  key=lambda r: r[at])
    if len(have) < 4: continue
    ys = [Y(r) for r in have]
    left, right = New(), adds(ys, New())
    for i, r in enumerate(have[:-1]):
      left  = add(left, ys[i])
      right = add(right, ys[i], w=-1)
      if r[at] == have[i+1][at]: continue
      n1, n2 = i + 1, len(have) - i - 1
      score = (n1*div(left) + n2*div(right)) / len(have)
      if best is None or score < best.score:
        cut  = (r[at] + have[i+1][at]) / 2
        best = o(at=at, cut=cut, score=score,
                 txt=f"{tbl.cols.names[at]} <= {cut:.3g}")
  return best

def tree(tbl, rows=None): # recursive supervised splitting
  Y, New, _ = yinfo(tbl)
  def go(rows):
    t = o(rows=rows, n=len(rows),
          ys=adds((Y(r) for r in rows), New()))
    if len(rows) > the.stop and (cut := split(tbl, rows)):
      left = lambda r: r[cut.at] != "?" and r[cut.at] <= cut.cut
      yes  = [r for r in rows if left(r)]
      no   = [r for r in rows if not left(r)]
      if 0 < len(yes) < len(rows):
        t.cut, t.yes, t.no = cut, go(yes), go(no)
    return t
  return go(tbl.rows if rows is None else rows)

def leaf(t, row): # walk row down to its leaf
  while hasattr(t, "cut"):
    c = t.cut
    ok = row[c.at] != "?" and row[c.at] <= c.cut
    t = t.yes if ok else t.no
  return t

def predict(tbl, t, row, k=None): # klass of k nearest in leaf
  rows = sorted(leaf(t, row).rows,
                key=lambda r: distx(tbl, row, r))[:k or the.k]
  at = tbl.cols.klass
  return mid(adds((r[at] for r in rows),
                  Col(tbl.cols.names[at])))


# ---------------------------------------------------------------
def test_disty():
  "sort rows by distance to goals; show best and worst"
  t = Tbl(csv(the.file))
  lst = sorted(t.rows, key=lambda r: disty(t, r))
  n = len(lst)
  assert disty(t, lst[0]) < disty(t, lst[-1])
  print(t.cols.names)
  for j in [0, 1, 2, n-3, n-2, n-1]:
    print(lst[j], "d %.2f" % disty(t, lst[j]))

def test_split():
  "show the split minimizing spread of the goals"
  t = Tbl(csv(the.file))
  cut = split(t, t.rows)
  print(cut.txt, "spread %.3f" % cut.score)

def test_tree():
  "print supervised tree with n and mean goal distance"
  t = Tbl(csv(the.file))
  t1 = tree(t)
  assert t1.n == sum(l.n for l in leaves(t1))
  show(t1, more=lambda n: " mu=%.2f" % mid(n.ys))
  print(len(leaves(t1)), "leaves")

def test_classify():
  "resubstitution accuracy on diabetes"
  t  = Tbl(csv("/Users/timm/gits/moot/classify/diabetes.csv"))
  t1 = tree(t)
  at = t.cols.klass
  acc = sum(predict(t, t1, r) == r[at]
            for r in t.rows) / len(t.rows)
  print("resub acc %.2f" % acc)
  assert acc > 0.7

eg = {
      "-disty": test_disty,
      "-split": test_split,
      "-tree": test_tree,
      "-classify": test_classify}

if __name__ == "__main__": main(meta(eg))
