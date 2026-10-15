# Punctured-E8 dual-basis decoders, p = 25..44: exhaustive exclusion (research notes)

Working scripts and notes, not a registered research package, and not independently audited. See NOTES.md for the
model, the exact class excluded, the argument, and the holes left open.

Scripts write intermediates (flats.npy, ftypes.npy, Cm.npy, minopts_root.pkl, adj44.pkl, viol.pkl, res2.pkl) into this
directory: setup.py (flats, types, label data) -> root_opts.py / minopts.py (kernels) -> family.py, dfs.py, dfs2.py
(exact family search; dfs2.py writes viol.pkl) -> analyze.py / analyze2.py (biorthogonality refinement at p = 36, 44).
d7*.py handle D7 kernels; sanity.py spot-checks the per-term rank inequality. Needs numpy (ortools only for ilp_root.py,
an early attempt that is not part of the argument).

Apache-2.0.
