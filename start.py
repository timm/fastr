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

def meta(eg): # give a demo table its --all and --egs entries
  def _all():
    "run all demos, in table order"
    for k, f in eg.items():
      if k[1] != "-": print("\n#", k); run(f)
  def _egs():
    "show the demo table"
    for k, f in eg.items():
      print("%-14s %s" % (k, (f.__doc__ or "").strip()))
  eg["--all"], eg["--egs"] = _all, _egs
  return eg

def main(eg): # run eg["-flag"]; else --setting gets a value
  for s, v in zip(sys.argv, sys.argv[1:] + [None]):
    if f := eg.get(s): run(f)
    elif s == "-h": print(__doc__)
    elif s.startswith("--"):
      k, eq, v1 = s[2:].partition("=")
      if hasattr(the, k): setattr(the, k, atom(v1 if eq else v))


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
eg = {
      "-the": test_the,
      "-atom": test_atom,
      "-csv": test_csv}

if __name__ == "__main__": main(meta(eg))
