# P112 two-star dictionary: decoder floor 64

**Result (EXACT, exhaustive, SUBCLASS).** The 16,556-entry symmetric two-star support dictionary of
[#5088](https://github.com/CrocSwap/integer-mult-bounds/pull/5088) contains **no masked P112 decoder of
cost ≤ 63**. Its minimum is **64**, attained by the coordinate decoder. #5088's own bounded MILP
(cost ≤ 63; lines, orthogonal pairs, frames) is therefore infeasible; their 90-second solve had timed out.
No new κ.

The proof needs only three facts, all from #5088: every dictionary price is ≥ 6; every orthogonal
eight-frame is met by ≥ 8 terms; every line is covered. Adding any number of images of the interleaved
atom43 term (either orientation) does not lower the floor either; that extension also uses the global
rank row (`Σ rank_Q ≥ 112`) in one case.

## Scope

Decoders `Σ_k (p_k q_k^T) ∘ Bd = I` (masked rank-1 model, DreamingOfClouds / PR717) whose terms have
symmetric support caps `S × S` with `S` in the #5088 dictionary, each occurrence priced at the cap
(`rank_F2 Bd[S,S]`), any complex coefficients, repeats allowed. Not covered: non-symmetric or
non-dictionary caps (the KFAM rank-1 lows, clique-like terms), other interleaved terms, coefficient
synthesis or any physical realization.

## Proof outline

See [PROOF.md](PROOF.md); each step is tagged PROVED or EXACT.

- **≤ 8 terms.** With 8 terms every term meets every frame; such supports all have price ≥ 8, so ≥ 64.
- **9 terms.** Each frame is missed at most once, so frame-missing terms are pairwise miss-disjoint.
  Rank-6 supports always share a missed frame, which forces ≥ 8 pairwise miss-disjoint price-7 supports.
  Up to symmetry such 8-sets exist only inside one orbit (orbit 11; 4,718,592 through its
  representative), and that orbit has no 9-set.
- **10 terms.** ≥ 7 distinct rank-6 terms with no frame missed three times: 1,248 such 7-sets up to
  symmetry. All 612,304 frame-feasible completions by 3 further supports leave a line uncovered.
- **Extension.** Frame-rank efficiency (dictionary ≤ 170 per unit cost, atom43 101/3) allows at most
  one atom; the remaining cases die by counting, a pair check and the global rank row.

## Reproduce

```sh
python3.12 -I -B research/p112-dictionary-floor/verify.py                  # about 2 min, 1 thread
python3.12 -I -B research/p112-dictionary-floor/verify.py --cross-check    # independent code, about 3 min
python3.12 -I -B research/p112-dictionary-floor/verify.py --full           # + the 4,718,592 count, about 2.5 min
```

Standard library only (Python ≥ 3.10). Expected marker: `PASS_P112_DICTIONARY_FLOOR_64`.
Saved outputs: [data/verify_primary.log](data/verify_primary.log),
[data/verify_crosscheck.log](data/verify_crosscheck.log), [data/verify_full.log](data/verify_full.log) and the JSON receipts beside them.

#5088's `verify.py` and `RESULT.json` are vendored unchanged under `vendor/masked-two-star-supports/`
(Apache-2.0) and checked against SHA256 pins at #5088 head `ce35592fc32351125435ecfd9f3f809bb898d312`.
To read them from Git instead, fetch the head and pass `--repo .`:

```sh
git fetch origin pull/5088/head && python3.12 -I -B research/p112-dictionary-floor/verify.py --repo .
```

`--pr-dir DIR` reads them from any directory. The verifier runs #5088's `geometry()`, `dictionary()` and
`cliques()` and requires its own rebuild to agree, including the dictionary SHA256 `9f62e399…dd50`.
`data/aut_generators.json` (11 permutations) is checked, not trusted: each must preserve orthogonality on
P112 and the dictionary. `data/atom43.json` gives atom43 as root vectors and levels.

## Caveats

- Exhaustive enumeration in Python, run twice with separate code; not a DRAT or formal-kernel certificate.
- `SOURCES.json` pins the vendored #5088 files to the unmerged #5088 head; the metadata check needs that
  commit fetched, and `--check-source-history` will fail until #5088 is merged.
- Support-level necessary conditions only. A cover at 64 says nothing about coefficients, and the floor
  says nothing outside the dictionary.

## Open next questions

1. **Cheap clique-like terms.**
   *Why it matters:* the frame argument depends on cost per frame met being ≤ 170. A cost-1 rank-1 term on
   a frame clique `F × F` meets 459 frames, and KFAM rank-1 lows and CT8-INT T-block Ferrers terms are
   similar, so the efficiency cap gives them no useful bound.
   *First step:* add the clique-supported terms of #5088's 262,928 cliques (priced at their cap) to the
   dictionary and rerun the 9/10-term enumeration with the per-term frame-rank in place of "meets".
   *Blocked by:* a price model for clique caps that matches the PR717 decoder cost; the term count grows
   fast and may need a new normalisation.
2. **Does the frame argument extend beyond the dictionary?**
   *Why it matters:* the proof used only price ≥ 6, frames ≥ 8 and lines ≥ 1, so it applies to any
   symmetric-cap family with the same price floor and frame census. A general version would be a floor
   for whole decoder families, not just one dictionary.
   *First step:* find the largest symmetric-cap family (for example all intersections of three stars)
   whose price ≥ 6 and whose pairwise shared-miss property for rank-6 caps still hold, and rerun `verify.py`
   on it.
   *Blocked by:* non-symmetric caps break the frame-meets count (rank on `F × F` then depends on both
   sides), and the 10-term completion grows with the family size.

## Credit

Dictionary, geometry, clique rows and the question (is there a cost-≤-63 cover?) are due to
**1337isnot1337 (Cinna Davis)**, #5088; the masked decoder model is DreamingOfClouds / PR717; the
frame-rank condition is stated in #5088 and its predecessors credited there. Their files are used
unchanged. Proof, enumeration and this package: Claude Opus 5.5 (Anthropic), via Claude Code, with an
independent audit of the enumeration whose code is folded into `--cross-check`. Apache-2.0.
