STATUS (Oct 14): window p=25..44 EXCLUDED for Model C, exactly 8 terms with biorthogonal normals (dual-basis decoders), exhaustive certificate. OPEN: bit model (P102's own model, 9 terms), m>8 terms, non-biorthogonal N.
# PUNCT2: exclusion of punctured-E8 decoders beating P102 (54/102) in p=25..44

## Class excluded
Model C of PUNCT/NOTES.md: P=R X, Q=R Y on the live lines K=L\S, sum_j x_j y_j^T = I_8 (8 terms, <x_i,y_j> = delta_ij c_j; N=XY^T diagonal), exact Gram on K,
cost = sum_j rank_F2 B[R_j&K, A_j&K], R_j = live lines not orth to x_j, A_j likewise for y_j. (Live K spans Q^8 for p<=44, so m>=8 automatically.)
Claim: for every S with |S| in 25..44, cost >= (54/102)(120-|S|) (no strict win; ties not excluded).

## Proof structure (all exact integer arithmetic / exhaustive)
1. Lemma (sanity.py, 3000 random checks, 0 violations; proof: Theorem 3 of NOTES): with Omega nondegenerate on F2^9, r_j = rank Omega(U_R,U_A) >= 9 - d(x_j) - d(y_j), d(x) = 9 - dim span labels(C_x & K).
2. d(x) <= d(f) for the rank-7 parabolic flat f containing x^perp (8 W-types, 440880 flats, setup.py). d(f)>=1 iff some functional z has O_z & C_f inside S. d(f) = log2(1+#active z).
3. Price of d>=2 kernels (cap 44; root_opts.py, d7c.py): types 2..7 (A7,E6A1,A1A6,D5A2,A4A3,A2A1A4): none <=44. D7: 64 options of exactly 42 lines per flat,
   no two D7 kernels within 2 lines of each other (69120 distinct), and a D7 kernel + 2 lines contains only 1 root d=2 option (d7b.py) => any S containing a D7 d=2 kernel has Top8 <= 10.
   Root flats: d=2 kernels 25..44 (10607 inclusion-minimal), d=3 kernels 37..44 (6573); no d=4.
4. Relaxation: cost >= 72 - D, D <= 2*Top8, Top8 <= 8 + E(S), E(S) = sum over root flats of (d-1)^+ (every non-root flat credited a free 1). No-win needs E <= 2 (p=25), 3 (26-29), 4 (30-33), 5 (34-36), 6 (37-40), 7 (41-44).
5. Exact max E over all S with |S|<=44 (family.py + dfs.py + dfs2.py): enumerate families of root flats with kernels, union<=44, fixing one flat to line 54 by W(E8) transitivity.
   Pair unions: two d=2 root kernels need union >=26 (orthogonal roots) / >=30 (non-orthogonal). Max E: 1 at 25, 2 for 26..35, **6 at 36**, 6 for 37..43, **8-12 at 44**.
   So the pure relaxation holds for every p except p=36 (E=6 > 5; 378 unions) and p=44 (E=8..12 > 7; 945 unions).
6. Refinement at those 1323 unions (analyze.py/analyze2.py): the 6 (resp. 12) d=2 roots R(S) computed exactly; sides with d>=2 must be roots in R; biorthogonality
   (<x_i,y_j>=0, i!=j; <x_j,y_j> != 0, distinct x's, distinct y's) over R plus free sides=1 gives max D = 22 (<28 needed) at p=36 and 24 (<32 needed) at p=44.
   (The 6 roots there cannot serve as both x and y sides: a root y_j orthogonal to all x_i, i!=j, forces y_j = x_j, and x=y terms cost more.)
=> no S in 25..44 admits cost < (54/102)(120-p). Together with PUNCT/NOTES (p<=24 bound 7/13 > .5294 via Theorem 4; p>=45 trivial? NOT shown here) the Model C 8-term dual-basis class is closed for p<=44.
Note: this excludes x < 54/102 (beating P102). Beating only P104 (x in [.5294,.5385)) is NOT covered here (thresholds would need recomputation with 56/104; E max 2 at 26..35 may violate them).

## Holes (OPEN, not claimed)
- Bit model (P102 itself is a 9-term bit-model decoder, R_j = Omega-orbits not hyperplane complements): not covered.
- m>8 (overcomplete) or non-diagonal N: biorthogonality was used only in step 6 (p=36, p=44); without it the relaxation permits D=28 (p=36) and 32 (p=44): relaxed cost 44 vs 44.47 and 40 vs 40.24.
- Free-side credit 1 for non-root flats and relaxed r>=9-dx-dy are upper bounds; true cost is >=; so nothing here constructs a decoder.
- p>=45 not treated (outside the window).
## Files: setup.py (flats), root_opts.py/minopts.py (kernels), family.py/dfs.py/dfs2.py (families), analyze*.py (biorthogonal refinement), d7*.py, sanity.py; results viol.pkl, res2.pkl.
