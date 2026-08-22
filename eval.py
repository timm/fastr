"""
eval: cross-val, confusion, effect size, significance, ranks.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license
"""
import sys; sys.dont_write_bytecode = True
from super import *
from bisect import bisect_right

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

def cliffs(xs, ys): # true if effect size is small
  gt = sum(x > y for x in xs for y in ys)
  lt = sum(x < y for x in xs for y in ys)
  return abs(gt - lt) / (len(xs)*len(ys) + 1/BIG) <= the.cliffs

def cohen(xs, ys): # true if means closer than cohen * pooled sd
  a, b = adds(xs, Num()), adds(ys, Num())
  sp = (((a[0]-1)*sd(a)**2 + (b[0]-1)*sd(b)**2) /
        (a[0] + b[0] - 2 + 1/BIG)) ** 0.5
  return abs(a[1] - b[1]) <= the.cohen * sp

def ks(xs, ys): # kolmogorov-smirnov accepts "same" at 95%?
  xs, ys = sorted(xs), sorted(ys)
  n, m = len(xs), len(ys)
  cdf = lambda a, v: bisect_right(a, v) / len(a)
  d = max(abs(cdf(xs, v) - cdf(ys, v)) for v in xs + ys)
  return d <= 1.36 * ((n + m)/(n*m)) ** 0.5

def same(xs, ys): # statistically indistinguishable?
  return cliffs(xs, ys) and ks(xs, ys)

def ranks(d): # {name:nums} --> {name:rank}; same() same rank
  out, rank, last = {}, 0, None
  order = sorted(d.items(),
                 key=lambda z: mid(adds(z[1], Num())))
  for name, xs in order:
    if last is not None and not same(last, xs): rank += 1
    out[name], last = rank, xs
  return out


# ---------------------------------------------------------------
def test_confuse():
  "confusion matrix stats on a small example"
  pairs = [("a","a")]*4 + [("b","a")] + \
          [("b","b")]*2 + [("a","b")]
  a = confuse(pairs)["a"]
  assert a.tp == 4 and a.fn == 1 and a.fp == 1 and a.tn == 2
  assert a.pd == 0.8 and abs(a.prec - 0.8) < 1e-9
  for c in confuse(pairs).values():
    print(c.label, "acc %.2f pd %.2f pf %.2f prec %.2f" %
          (c.acc, c.pd, c.pf, c.prec))

def test_stats():
  "same() accepts tiny shifts, rejects big ones"
  a = [rand() for _ in range(40)]
  b = [x + 0.05 for x in a]
  c = [x + 2    for x in a]
  assert same(a, b) and not same(a, c)
  assert cohen(a, b) and not cohen(a, c)
  print("a=b:", same(a, b), "a=c:", same(a, c))

def test_ranks():
  "same rank for same distributions, else new rank"
  d = {"x1": [rand()       for _ in range(30)],
       "x2": [rand()       for _ in range(30)],
       "x3": [rand() + 0.5 for _ in range(30)],
       "x4": [rand() + 2   for _ in range(30)]}
  r = ranks(d)
  print(r)
  assert r["x1"] == r["x2"] < r["x3"] < r["x4"]

def test_xval():
  "cross-val classification on diabetes and soybean"
  for f in ["diabetes", "soybean"]:
    print("\n##", f)
    t = Tbl(csv(f"/Users/timm/gits/moot/classify/{f}.csv"))
    k, pairs = t.cols.klass, []
    for train, test in xval(t.rows, n=3, m=5):
      tr = clone(t, train)
      t1 = tree(tr)
      for row in test:
        pairs += [(predict(tr, t1, row), row[k])]
    for c in sorted(confuse(pairs).values(),
                    key=lambda c: -c.pd):
      print("%25s acc %.2f pd %.2f pf %.2f prec %.2f" %
            (c.label, c.acc, c.pd, c.pf, c.prec))

def test_opt():
  "optimizer: tree picks beat random picks"
  t = Tbl(csv(the.file))
  best, rnd = Num(), Num()
  for train, test in xval(t.rows, n=3, m=10):
    tr = clone(t, train)
    t1 = tree(tr)
    score = lambda r: mid(leaf(t1, r).ys)
    picks = sorted(test, key=score)[:the.check]
    some  = sample(test, the.check)
    best  = add(best, min(disty(tr, r) for r in picks))
    rnd   = add(rnd,  min(disty(tr, r) for r in some))
  lo = min(disty(t, r) for r in t.rows)
  print("true best %.2f | best of %s picked %.2f | random %.2f" %
        (lo, the.check, mid(best), mid(rnd)))

eg = {
      "-confuse": test_confuse,
      "-stats": test_stats,
      "-ranks": test_ranks,
      "-xval": test_xval,
      "-opt": test_opt}

if __name__ == "__main__": main(meta(eg))
