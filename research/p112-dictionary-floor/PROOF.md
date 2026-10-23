# Proof: the #5088 dictionary has no masked P112 decoder of cost at most 63

Tags: **PROVED** = short argument written here; **EXACT** = finite computation that `verify.py`
reruns (default mode) and reruns again with separate code (`--cross-check`). Nothing is sampled.
The computations are exhaustive enumerations in Python, not DRAT or formal-kernel certificates.

## Scope (SUBCLASS)

The decoder model is the masked rank-1 model of DreamingOfClouds / PR717, as used in #5088.
Let `Bd` be the 112 x 112 0/1 matrix on the P112 lines with `Bd[a,b] = 1` iff `a = b` or `a ⟂ b`.
A decoder is a finite multiset of terms `T_k = (p_k q_k^T) ∘ Bd` with `Σ_k T_k = I` (complex
coefficients, cancellation allowed). The covered class is:

- every term has a **symmetric support cap** `S_k × S_k`, where `S_k` is one of the 16,556 supports of
  the #5088 two-star dictionary (`supp p_k, supp q_k ⊆ S_k`);
- each occurrence of a term costs the cap's price `price(S) = rank_F2 Bd[S,S]` (the #5088 price);
- multiplicities are allowed and each occurrence pays.

Not covered: non-symmetric caps, caps outside the dictionary (in particular the rank-1 KFAM low terms
and clique-like terms), interleaved/Ferrers terms other than atom43 (section 6), and anything about
coefficients, physical COPY realization or κ.

Notation: a **frame** is one of the 1,260 orthogonal 8-cliques of P112; a support **meets** a frame `F`
if `S ∩ F ≠ ∅`, otherwise it **misses** `F`. `m` is the number of terms (with multiplicity) and
`a_p` the number of terms of price `p`.

## 1. The three necessary conditions (PROVED; due to #5088)

1. **Price ≥ 6.** Every dictionary support has price 6, 7, 8 or 9 (EXACT: 3,780 / 3,360 / 6,392 / 3,024).
2. **Frames ≥ 8.** On `F × F` every entry of `Bd` is 1, so `T_k` restricted to `F × F` is
   `p_k|_F q_k|_F^T`: rank ≤ 1, and 0 when `S_k` misses `F`. The restriction of `I` has rank 8, so by
   subadditivity at least 8 terms (with multiplicity) meet `F`.
3. **Lines ≥ 1.** `I[u,u] = 1` needs some term with `u ∈ S_k`.

#5088's MILP used rows 2 and 3 plus orthogonal pairs ≥ 2; the argument below needs only 1-3, so it
also shows that their 4,676-row MILP with cost ≤ 63 is infeasible (their 90 s run timed out).

## 2. Data (EXACT)

- The dictionary is rebuilt by #5088's own `dictionary()` (bytes pinned to head `ce35592f`) and by
  separate code; both give the same 16,556 masks and prices, SHA256 `9f62e399…dd50` = #5088 `RESULT.json`.
  The price/frames-met census equals #5088's `atom_frame_census`.
- `data/aut_generators.json` holds 11 permutations of the root lines. `verify.py` checks that each
  fixes the 8 removed lines and preserves orthogonality on P112; hence each preserves `Bd`, the frames,
  the dictionary (checked) and prices. Only "these generate a group of automorphisms" is used.
  The supports fall into 16 orbits (numbered by smallest member in sorted-mask order):

  | orbit | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | size | 224 | 448 | 224 | 2688 | 5376 | 2688 | 2688 | 448 | 448 | 28 | 168 | 448 | 448 | 112 | 112 | 8 |
  | price | 7 | 8 | 9 | 7 | 8 | 9 | 6 | 6 | 6 | 6 | 6 | 7 | 8 | 8 | 9 | 8 |

- **Normalisation (PROVED).** Conditions 1-3 are invariant under the group, so any decoder can be
  moved so that a chosen sub-multiset contains the representative of its smallest orbit and all its
  other members lie in orbits ≥ that one.

## 3. At most 8 terms (PROVED + EXACT)

`m < 8` violates condition 2. If `m = 8` every term meets every frame. EXACT: the 3,480 supports that
meet all 1,260 frames all have price ≥ 8, so the cost is ≥ 64. If `m ≥ 11` the cost is ≥ 66.
So a cost-≤-63 decoder has `m = 9` or `m = 10`.

Control (EXACT): the coordinate decoder (the 8 supports of orbit 15) costs 64, meets every frame
exactly 8 times and covers every line; dropping one of its terms breaks a frame.

## 4. Nine terms (PROVED + EXACT)

With `m = 9` each frame is missed by at most one term (counting multiplicity). So terms that miss a
frame have multiplicity 1 and pairwise disjoint miss sets.

- EXACT: every price-6/7 support misses some frame; all 7,142,310 pairs of rank-6 supports share a
  missed frame (exhaustive). Hence `a_6 ≤ 1`.
- Cost ≤ 63 gives `Σ(p − 7) = −a_6 + a_8 + 2a_9 ≤ 0`, so `a_8 + 2a_9 ≤ a_6 ≤ 1`: the profiles are
  `(a_6,a_7,a_8,a_9) = (0,9,0,0), (1,8,0,0), (1,7,1,0)`. Each contains ≥ 8 pairwise miss-disjoint
  price-6/7 supports, at most one of price 6.
- EXACT: normalise such an 8-set (section 2). For every price-6/7 orbit `j ≠ 11` the largest
  pairwise miss-disjoint set containing `rep_j`, drawn from orbits ≥ `j` with at most one rank-6 member,
  has size ≤ 7. So the 8-set lies inside orbit 11 (price 7), hence `a_6 = 0`, the profile is
  `(0,9,0,0)` and all 9 terms are pairwise miss-disjoint members of orbit 11.
- EXACT: the largest pairwise miss-disjoint subset of orbit 11 containing `rep_11` has size 8.
  No 9-set exists: `m = 9` is impossible.

(Not needed for the proof: orbit 11 has 4,718,592 such 8-sets through `rep_11`; `--full` recounts them.)

## 5. Ten terms (PROVED + EXACT)

With `m = 10` each frame is missed at most twice. Cost ≤ 63 gives `Σ(p − 6) ≤ 3`, so `a_6 ≥ 7`.
The rank-6 terms are distinct: a repeated rank-6 support together with any other rank-6 term
(section 4: they share a missed frame) misses some frame three times.

- EXACT: up to the group (smallest orbit's representative included, the rest from orbits ≥ it), the
  sets of 7 distinct rank-6 supports with no frame missed three times number **1,248**, all with
  smallest orbit 7.
- EXACT: complete each with 3 arbitrary dictionary supports (any orbit, repeats allowed, total price
  ≤ 63 − 42 = 21, which covers all seven price profiles) such that no frame is missed three times:
  612,304 completions, and **every one leaves a line uncovered**. `m = 10` is impossible.

**Theorem (EXACT, SUBCLASS).** In the class defined under Scope, every decoder costs at least 64, and the
coordinate decoder attains 64.

## 6. Extension: dictionary + atom43 (PROVED + EXACT)

atom43 is the interleaved (Ferrers) term `(q p^T) ∘ S` on 4 target × 5 source lines (`data/atom43.json`,
`S[r][c] = [level(c) ≤ level(r)]`), cost `rank_F2(Bd ∘ S) = 3`, `rank_Q(Bd ∘ S) = 4`. Any image under
the group, in either orientation, is allowed, any number of times.

- **Frame-rank lemma (PROVED).** For a term `T`, let `fr(T) = Σ_F rank T|_{F×F}`. By condition 2's
  argument, `Σ_k fr(T_k) ≥ 8 · 1,260 = 10,080`. A dictionary term has `fr ≤` (frames met);
  atom43 has `fr ≤ Σ_F rank S[F∩R, F∩C]`.
- EXACT: dictionary frames-met / price ≤ 170 (orbits 7, 9); atom43 has rank 1 on exactly 101 frames and
  0 on the rest (`fr ≤ 101`; transposing changes nothing). With `k` atom images and dictionary cost
  `c_D ≤ 63 − 3k`: `170(63 − 3k) + 101k < 10,080` for `k ≥ 2`. So `k ≤ 1`; `k = 0` is the theorem.
- `k = 1` (`c_D ≤ 60`, atom fixed by normalisation; `A` = its 101 frames). A frame is missed by at
  most `m_D − 8 + [F ∈ A]` dictionary terms.
  - `m_D ≤ 7`: a frame outside `A` needs 8 dictionary hits. `m_D ≥ 11`: `c_D ≥ 66`.
  - `m_D = 8`: dictionary misses occur only on `A`, at most once each (≤ 101 in total). EXACT:
    price-≤7 supports miss ≥ 90 frames, so at most one term has price ≤ 7 and `c_D ≥ 6 + 7·8 = 62`.
  - `m_D = 9`: `c_D ≤ 60` gives `a_6 ≥ 3`. Two rank-6 terms may share a missed frame only inside `A`.
    EXACT: no pair of rank-6 supports (repeats included) has its shared misses inside `A`. Dead.
  - `m_D = 10`: all ten terms have price 6. **Rank row (PROVED, a fourth ingredient):** `T_k` has rank
    ≤ `rank_Q Bd[S_k,S_k]` (`= 7` for every price-6 support, EXACT) and the atom ≤ 4, while `I` has rank
    112: `10·7 + 4 = 74 < 112`. Dead.

**Corollary (EXACT, SUBCLASS).** Dictionary terms plus any number of atom43 images still cost ≥ 64.
The extension uses the global rank row in the case `m_D = 10`; the main theorem does not.

## Cross-check (`--cross-check`)

Separate code for every EXACT step: prices by explicit row-span enumeration, frames by a
pairwise-orthogonality recursion, miss sets by direct scan, orbits by union-find with a different
orbit labelling (sorted by price and size), the rank-6 pair check by a direct double loop, maximum
miss-disjoint sets by branch and bound with a colouring bound, triple-free 7-sets in the other
labelling, and the 10-term completion by a line-driven search at budget 21 (control: budget 45 finds a
covering completion, so the line test is not vacuous). This mode adapts the independent auditor's code.
