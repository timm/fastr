# fastr

Very simple cluster, contrast, classify, regress, multi-objective.

Each file builds on the last. Run `python3 FILE.py -all` for that file's demos,
`-h` for help, `--flag value` to change settings (see start.py).

- start.py    : config, atoms, csv rows, test-runner CLI
- tbl.py      : Num, Sym, Cols, Tbl (incremental column summaries)
- contrast.py : ranges that most separate two row sets (the.bins cuts of -3sd..+3sd)
- unsuper.py  : distx + fastmap recursive bi-clustering; contrast names the splits
- super.py    : disty + supervised splits (minimize spread of goals); classify, regress
- eval.py     : xval, confusion, cliffs delta, ks test, cohen d, same-rank ranking
