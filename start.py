"""
start: config, atoms, csv rows, and a tiny test-runner CLI.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license

Options:

  -h             show help
  --p=2          minkowski coefficient
  --stop=4       stopping rule for recursive tree generation
  --few=256      sub-sample size for pole finding
  --k=5          nearest neighbors used in a leaf
  --check=5      optimization: how many top picks to evaluate
  --bins=10      number of bins for discretization
  --cliffs=.197  cliffs delta: max effect size for "same"
  --cohen=.35    cohen d: max mean separation for "same"
  --seed=1234    random number generation
  --file=/Users/timm/gits/moot/optimize/misc/auto93.csv
"""
import re, sys; sys.dont_write_bytecode = True
from math import exp, log
from types import SimpleNamespace as o
from random import choice, sample, seed, shuffle, random as rand

BIG = 1e32

def atom(s): # string --> atom
  try: return int(s)
  except ValueError:
    try: return float(s)
    except ValueError: return s.strip()

the = o(**{k: atom(v)
           for k, v in re.findall(r"(\w+)=(\S+)", __doc__)})

def csv(f): # csv file strings --> rows of atoms
  with open(f, encoding="utf-8") as fs:
    for s in fs:
      if s := s.strip(): yield [atom(x) for x in s.split(",")]


# ---------------------------------------------------------------
def run(f, *v):
  try: seed(the.seed); f(*v)
  except Exception: import traceback; traceback.print_exc()

def demos(g): # a demo = local test_ function with a docstring
  return {k[5:]: f for k, f in list(g.items())
          if k.startswith("test_") and f.__doc__
          and f.__module__ == g["__name__"]}

def main(g): # single dash = local demo, double dash = setting
  for s, v in zip(sys.argv, sys.argv[1:] + [None]):
    if s == "-h":
      print(__doc__ + "\nDemos:\n")
      for k, f in demos(g).items():
        print("  -%-13s %s" % (k, f.__doc__.strip()))
    elif s == "-all":
      for k, f in demos(g).items(): print("\n#", k); run(f)
    elif s.startswith("--"): # --key=value or --key value
      k, eq, v1 = s[2:].partition("=")
      if hasattr(the, k): setattr(the, k, atom(v1 if eq else v))
    elif f := demos(g).get(s[1:].replace("-", "_")):
      run(f)


# ---------------------------------------------------------------
def test_the():
  "show current settings"
  print(the)

def test_atom():
  "strings coerce to ints, floats, or stripped strings"
  assert atom("2") == 2 and atom("2.1") == 2.1
  assert atom(" a ") == "a"
  print("'2' ->", atom("2"), "| '2.1' ->", atom("2.1"),
        "| ' a ' ->", repr(atom(" a ")))

def test_csv():
  "csv reader finds many rows in the.file"
  rows = list(csv(the.file))
  assert len(rows) > 100
  print(len(rows), "rows; first:", rows[0])

seed(the.seed)
if __name__ == "__main__": main(globals())
