"""
contrast: ranges that most separate two groups of rows.
X ranges come from the parent table; one XY per column counts
which group lands in each range. Contrast is just counting.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license
"""
import sys; sys.dont_write_bytecode = True
from tbl import *

def score(z): return z.score

def discretize(d0, d1, ordered=False): # two count dicts
  "Yield spans scored b^2/(b+r); b,r = frequency in d0,d1."
  n0 = sum(d0.values()) + 1/BIG
  n1 = sum(d1.values()) + 1/BIG
  ks = sorted(set(d0) | set(d1))
  B, R, m = [0], [0], len(ks)
  for k in ks:
    B += [B[-1] + d0.get(k, 0)/n0]
    R += [R[-1] + d1.get(k, 0)/n1]
  for i in range(m):
    for j in (range(i + 1, m + 1) if ordered else [i + 1]):
      if (i, j) != (0, m): # whole column selects nothing
        yield _span(ks[i], ks[j-1], B[j]-B[i], R[j]-R[i],
                    i == 0, j == m)

def _span(lo, hi, b, r, first, last): # the scorer: b^2/(b+r)
  return o(lo=lo, hi=hi, first=first, last=last,
           score=b*b / (b + r + 1/BIG))

def contrasts(tbl, rows0, rows1): # best range splitting groups
  cs = ranges(tbl)
  def go(at):
    xy = XY(cs.get(at))
    for y, rows in ((0, rows0), (1, rows1)):
      for r in rows: add(xy, (r[at], y))
    d0 = {b: d.get(0, 0) for b, d in xy.seen.items()}
    d1 = {b: d.get(1, 0) for b, d in xy.seen.items()}
    z = max((z for d in (discretize(d0, d1, at in cs),
                         discretize(d1, d0, at in cs))
             for z in d), key=score, default=None)
    if z:
      z.at, name = at, tbl.cols.names[at]
      if at not in cs:
        z.txt = f"{name} == {z.lo}"
      else:
        z.lo = -BIG if z.first else cs[at][z.lo]
        z.hi =  BIG if z.last  else cs[at][z.hi + 1]
        z.txt = (f"{name} <= {z.hi:.3g}" if z.lo == -BIG else
                 f"{name} >= {z.lo:.3g}" if z.hi == BIG else
                 f"{name} in {z.lo:.3g}..{z.hi:.3g}")
    return z
  return max((z for at in tbl.cols.x if (z := go(at))),
             key=score, default=None)

def selects(cut, row): # does row fall inside a contrast range?
  v = row[cut.at]
  return v != "?" and cut.lo <= v <= cut.hi


# ---------------------------------------------------------------
def test_discretize():
  "two bin-count dicts; best span covers the first"
  d0 = {0: 30, 1: 50, 2: 20}
  d1 = {6: 40, 7: 40, 8: 20}
  z = max(discretize(d0, d1, ordered=True), key=score)
  assert z.score > 0.9 and z.lo == 0 and z.hi == 2
  print("d0 in bins 0-2, d1 in bins 6-8:",
        "span %s..%s score %.2f" % (z.lo, z.hi, z.score))

def test_contrasts():
  "name the 50:50 halves of fastmap's top-level split"
  import unsuper
  t = Tbl(csv(the.file))
  x = lambda r1, r2: unsuper.distx(t, r1, r2)
  _, rows = unsuper.fastmap(t.rows[:], x)
  n = len(rows) // 2
  r0, r1 = rows[:n], rows[n:]
  z = contrasts(t, r0, r1)
  p0 = sum(selects(z, r) for r in r0) / len(r0)
  p1 = sum(selects(z, r) for r in r1) / len(r1)
  print("best: %s; selects %.0f%% vs %.0f%%" %
        (z.txt, 100*p0, 100*p1))
  assert abs(p0 - p1) > 0.2

def test_halve():
  "one line: best range splitting fastmap's two halves"
  import unsuper
  t = Tbl(csv(the.file))
  x = lambda r1, r2: unsuper.distx(t, r1, r2)
  _, rows = unsuper.fastmap(t.rows[:], x)
  n = len(rows) // 2
  z = contrasts(t, rows[:n], rows[n:])
  print("%-45s %.2f %s" % (the.file.split("/")[-1],
                           z.score if z else 0,
                           z.txt if z else "no split"))

if __name__ == "__main__": main(globals())
