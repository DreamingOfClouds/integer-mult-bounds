#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact floor 64 for masked P112 decoders built from the #5088 two-star support dictionary.

Claim (SUBCLASS, EXACT, exhaustive): no multiset of rank-1 masked terms whose symmetric support caps
come from the 16,556-entry #5088 dictionary, each priced at its cap, has total price <= 63 while
(i) every one of the 1,260 orthogonal eight-frames is met by >= 8 terms and (ii) every line is covered.
The coordinate decoder (cost 64) meets both, so the dictionary minimum is 64.  Extension: adding any
number of images of the interleaved atom43 term (either orientation) does not help.

Standard library only (Python >= 3.10).  Modes:
  default          rebuild + all proof steps (about 3 minutes, one thread)
  --cross-check    independent second implementation of every step (different algorithms)
  --full           also count the 8-sets of pairwise miss-disjoint orbit-11 supports (not load-bearing)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent
PR_HEAD = "ce35592fc32351125435ecfd9f3f809bb898d312"          # CrocSwap/integer-mult-bounds PR #5088 head
PR_DIR = "research/masked-two-star-supports"
PINS = {"verify.py": "bc7e41ed70e095dbabd19ab4f8836bd299a64ee12872b7debc3eb33f04df7228",
        "RESULT.json": "d650e675ea6457f8fb165cc198db72f3f2cfe358271c3f020f6440d45792e9a0"}
DICT_SHA = "9f62e399532aaff3be299a6efc23450d7af04cd2351fb56ff91d2f5714c9dd50"
BUDGET = 63
TRIPLES = [t for t in product((6, 7, 8, 9), repeat=3) if t[0] <= t[1] <= t[2] and sum(t) <= BUDGET - 42]
T0 = time.time()
LOG = []


def say(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    print(line, flush=True)
    LOG.append(line)


def need(ok, msg):
    if not ok:
        raise AssertionError(msg)


def bits(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x ^= b


# ----------------------------------------------------------------------------------------------
# #5088 sources (pinned bytes)
# ----------------------------------------------------------------------------------------------
def read_pr_file(name, pr_dir, repo):
    if pr_dir is not None:
        data, where = (Path(pr_dir) / name).read_bytes(), f"--pr-dir/{name}"
    elif repo is not None:
        data = subprocess.run(["git", "-C", str(repo), "show", f"{PR_HEAD}:{PR_DIR}/{name}"],
                              check=True, capture_output=True).stdout
        where = f"git {PR_HEAD[:12]}:{PR_DIR}/{name}"
    else:
        data, where = (HERE / "vendor" / "masked-two-star-supports" / name).read_bytes(), f"vendor/{name}"
    h = hashlib.sha256(data).hexdigest()
    need(h == PINS[name], f"{where}: sha256 {h} != pinned {PINS[name]}")
    return data, where


def load_pr5088(pr_dir, repo):
    src, where = read_pr_file("verify.py", pr_dir, repo)
    res, _ = read_pr_file("RESULT.json", pr_dir, repo)
    mod = {"__name__": "pr5088_verify", "__file__": str(HERE / "vendor" / "masked-two-star-supports" / "verify.py")}
    exec(compile(src, f"pr5088:{PR_DIR}/verify.py", "exec"), mod)
    say(f"#5088 sources: {where} (+ RESULT.json), sha256 pinned to head {PR_HEAD}")
    return mod, json.loads(res)


# ----------------------------------------------------------------------------------------------
# geometry, dictionary, frames (own code; compared with #5088's)
# ----------------------------------------------------------------------------------------------
def own_geometry():
    full = set()
    for i, j in combinations(range(8), 2):
        for a, b in product((-2, 2), repeat=2):
            r = [0] * 8; r[i], r[j] = a, b; full.add(tuple(r))
    for v in product((-1, 1), repeat=8):
        if v.count(-1) % 2 == 0:
            full.add(v)
    need(len(full) == 240, "E8 root count")
    lines = sorted(r for r in full if next(x for x in r if x) > 0)        # one representative per line
    dot = [[sum(a * b for a, b in zip(x, y)) for y in lines] for x in lines]
    frame = [lines.index(tuple(2 if k == i else s * 2 if k == i + 1 else 0 for k in range(8)))
             for i in range(0, 8, 2) for s in (-1, 1)]
    live = [i for i in range(120) if i not in frame]
    LIVE = sum(1 << i for i in live)
    orth = [sum(1 << j for j in live if dot[i][j] == 0) for i in range(120)]   # orthogonal live lines
    bd = [orth[i] | (1 << i) for i in range(120)]                            # Bd row: same line or orthogonal
    return dict(roots=lines, dot=dot, frame_roots=frame, live=live, LIVE=LIVE, orth=orth, bd=bd)


def f2rank_low(rows):
    """F2 rank, lowest-bit pivots (the #5088 verifier uses highest-bit pivots)."""
    basis = {}
    for r in rows:
        while r:
            p = r & -r
            if p in basis:
                r ^= basis[p]
            else:
                basis[p] = r
                break
    return len(basis)


def price_of(G, s):
    return f2rank_low([G["bd"][i] & s for i in bits(s)])


def own_dictionary(G):
    LIVE, dot = G["LIVE"], G["dot"]
    stars = [sum(1 << j for j in range(120) if dot[i][j] != 0) & LIVE for i in range(120)]
    base = set(stars) | {s & ~(1 << i) for i, s in enumerate(stars)}
    sup = base | {a & b for a, b in combinations(base, 2)}
    sup.discard(0)
    sup = sorted(sup)
    return sup, [price_of(G, s) for s in sup], len(base)


def dict_hash(sup, pr):
    return hashlib.sha256("".join(f"{s:x}\t{r}\n" for s, r in zip(sup, pr)).encode()).hexdigest()


def frames_bitset(G):
    out = []
    def rec(ch, cand, k):
        if k == 8:
            out.append(ch); return
        while cand:
            b = cand & -cand; i = b.bit_length() - 1; cand ^= b
            rec(ch | b, cand & G["orth"][i], k + 1)
    rec(0, G["LIVE"], 0)
    return out


def qrank(M):
    """exact rank over Q (fraction-free Bareiss)."""
    A = [list(map(int, r)) for r in M]
    if not A or not A[0]:
        return 0
    n, m = len(A), len(A[0]); r = 0; prev = 1
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(r + 1, n):
            A[i] = [(A[r][c] * A[i][j] - A[i][c] * A[r][j]) // prev for j in range(m)]
        prev = A[r][c]; r += 1
        if r == n:
            break
    return r


def bd_block(G, rows, cols):
    return [[1 if (G["bd"][a] >> b) & 1 else 0 for b in cols] for a in rows]


# ----------------------------------------------------------------------------------------------
# automorphisms and orbits
# ----------------------------------------------------------------------------------------------
def permute(mask, g):
    r = 0
    for i in bits(mask):
        r |= 1 << g[i]
    return r


def check_generators(G, gens):
    live = set(G["live"])
    for g in gens:
        need(sorted(g) == list(range(120)), "generator is a permutation")
        need(all(g[i] == i for i in G["frame_roots"]), "generator fixes the removed frame lines")
        need(all(g[i] in live for i in live), "generator preserves the live lines")
        need(all(permute(G["orth"][i], g) == G["orth"][g[i]] for i in G["live"]), "generator preserves orthogonality")


def orbits_bfs(sup, gens):
    idx = {s: k for k, s in enumerate(sup)}
    orb = [-1] * len(sup); reps = []
    for k in range(len(sup)):
        if orb[k] >= 0:
            continue
        o = len(reps); reps.append(k); orb[k] = o; stack = [k]
        while stack:
            a = stack.pop()
            for g in gens:
                im = idx.get(permute(sup[a], g))
                need(im is not None, "dictionary is invariant under the generators")
                if orb[im] < 0:
                    orb[im] = o; stack.append(im)
    return orb, reps


# expected orbit census (orbit = numbered by smallest member in sorted-mask order): size, |S|, price
ORBITS = [(224, 24, 7), (448, 25, 8), (224, 26, 9), (2688, 26, 7), (5376, 27, 8), (2688, 28, 9), (2688, 22, 6),
          (448, 24, 6), (448, 24, 6), (28, 24, 6), (168, 20, 6), (448, 28, 7), (448, 29, 8), (112, 52, 8),
          (112, 53, 9), (8, 56, 8)]


# ----------------------------------------------------------------------------------------------
# primary proof
# ----------------------------------------------------------------------------------------------
def primary(args):
    R = {}
    mod, res5088 = load_pr5088(args.pr_dir, args.repo)
    roots, dots, live5088, pairing, adj5088 = mod["geometry"]()
    atoms = mod["dictionary"](dots, live5088, pairing)
    G = own_geometry()
    need(roots == G["roots"] and live5088 == G["LIVE"], "own geometry == #5088 geometry")
    need(all(adj5088[i] == G["orth"][i] for i in range(120)), "own orthogonality == #5088 adjacency")
    sup, pr, nbase = own_dictionary(G)
    need(nbase == 232 and len(sup) == 16556, "232 base supports, 16,556 dictionary entries")
    need([(s, r) for s, r in atoms] == list(zip(sup, pr)), "own dictionary and prices == #5088 dictionary()")
    h = dict_hash(sup, pr)
    need(h == DICT_SHA == res5088["dictionary_sha256"], "dictionary sha256")
    hist = {str(p): pr.count(p) for p in sorted(set(pr))}
    need(hist == res5088["rank_histogram"], "rank histogram == #5088 RESULT.json")
    say(f"dictionary: 16,556 supports, sha256 {h[:16]}... == #5088 RESULT.json; prices {hist}")
    K = len(sup); LIVE = G["LIVE"]

    FR = frames_bitset(G)
    theirs = sorted(c for c, s in mod["cliques"](adj5088, live5088, {8}) if s == 8)
    need(len(FR) == 1260 and sorted(FR) == theirs, "1,260 frames, == #5088 cliques()")
    for f in FR:                                   # frames are all-ones blocks of Bd
        need(all((G["bd"][i] & f) == f for i in bits(f)), "frame block of Bd is all ones")
    NF = len(FR); FULLF = (1 << NF) - 1
    FI = {u: 0 for u in G["live"]}
    for j, f in enumerate(FR):
        for u in bits(f):
            FI[u] |= 1 << j
    hit = []
    for s in sup:
        x = 0
        for u in bits(s):
            x |= FI[u]
        hit.append(x)
    miss = [FULLF ^ x for x in hit]
    census = {}
    for p, x in zip(pr, hit):
        key = f"{p}:{x.bit_count()}"; census[key] = census.get(key, 0) + 1
    need(dict(sorted(census.items())) == res5088["atom_frame_census"], "price:frames-met census == #5088")
    say(f"frames: 1,260 (== #5088); price:frames-met census == #5088 RESULT.json")

    gens = json.loads((HERE / "data" / "aut_generators.json").read_text())["generators"]
    check_generators(G, gens)
    for g in gens:
        fs = set(FR); need(all(permute(f, g) in fs for f in FR), "frames invariant")
    orb, reps = orbits_bfs(sup, gens)
    table = [(orb.count(o), sup[k].bit_count(), pr[k]) for o, k in enumerate(reps)]
    need(table == ORBITS, f"orbit census {table}")
    say(f"group: {len(gens)} checked automorphism generators; 16 support orbits " +
        " ".join(f"{o}:{n}/p{p}" for o, (n, _, p) in enumerate(table)))

    # ---- step 1: m <= 8 -------------------------------------------------------------------
    need(min(pr) == 6, "every dictionary price >= 6")
    full_hit = [k for k in range(K) if miss[k] == 0]
    need(min(pr[k] for k in full_hit) == 8, "supports meeting all 1,260 frames have price >= 8")
    need(all(miss[k] for k in range(K) if pr[k] <= 7), "every price-6/7 support misses a frame")
    say(f"step 1 (m<=8): min price 6 => m <= 10; {len(full_hit)} supports meet every frame, all price >= 8 => m = 8 costs >= 64")

    # positive control: coordinate decoder
    coord = [sum(1 << j for j in G["live"] if G["dot"][f][j] != 0) for f in G["frame_roots"]]
    idx = {s: k for k, s in enumerate(sup)}
    ck = [idx[c] for c in coord]
    need(sum(pr[k] for k in ck) == 64, "coordinate decoder costs 64")
    need(all(sum(1 for k in ck if not miss[k] >> j & 1) >= 8 for j in range(NF)), "coordinate decoder meets frames 8x")
    cov = 0
    for k in ck: cov |= sup[k]
    need(cov == LIVE, "coordinate decoder covers every line")
    need(any(sum(1 for k in ck[:-1] if not miss[k] >> j & 1) < 8 for j in range(NF)), "7 coordinates fail a frame")
    say("control: coordinate decoder (8 x orbit 15) costs 64 and meets every frame 8 times; dropping one fails")

    # ---- step 2: m = 9 ------------------------------------------------------------------------
    R6 = [k for k in range(K) if pr[k] == 6]
    HR6 = [0] * NF
    for a, k in enumerate(R6):
        for j in bits(hit[k]):
            HR6[j] |= 1 << a
    for a, k in enumerate(R6):
        x = (1 << len(R6)) - 1
        for j in bits(miss[k]):
            x &= HR6[j]
            if not x:
                break
        need(x == 0, "two rank-6 supports with disjoint miss sets")
    say(f"step 2a: all {len(R6) * (len(R6) - 1) // 2:,} pairs of rank-6 supports share a missed frame (exhaustive)")
    V = [k for k in range(K) if pr[k] in (6, 7)]; nV = len(V)
    HV = [0] * NF
    for a, k in enumerate(V):
        for j in bits(hit[k]):
            HV[j] |= 1 << a
    adj = []
    for a, k in enumerate(V):
        x = ((1 << nV) - 1) ^ (1 << a)
        for j in bits(miss[k]):
            x &= HV[j]
        adj.append(x)
    is6 = sum(1 << a for a, k in enumerate(V) if pr[k] == 6)
    vpos = {k: a for a, k in enumerate(V)}

    def max_clique(r, P, cap):
        """largest set of pairwise miss-disjoint supports containing r (<= one rank-6), capped at cap."""
        best = [1]
        def go(size, P, n6):
            if size > best[0]:
                best[0] = size
            if best[0] >= cap:
                return
            while P:
                if size + P.bit_count() <= best[0]:
                    return
                b = P & -P; i = b.bit_length() - 1; P ^= b
                k6 = (is6 >> i) & 1
                if n6 + k6 > 1:
                    continue
                go(size + 1, P & adj[i], n6 + k6)
        go(1, P, (is6 >> r) & 1)
        return best[0]

    m9 = {}
    for j, rk in enumerate(reps):
        if pr[rk] not in (6, 7):
            continue
        r = vpos[rk]
        P = adj[r] & sum(1 << a for a, k in enumerate(V) if orb[k] >= j)
        if j == 11:
            m9[j] = max_clique(r, P, 9)
            need(m9[j] == 8, "orbit 11: a 9-set of pairwise miss-disjoint supports exists")
        else:
            m9[j] = max_clique(r, P, 8)
            need(m9[j] < 8, f"orbit {j}: an 8-set of pairwise miss-disjoint supports exists")
    say(f"step 2b (m=9): max pairwise miss-disjoint set (<=1 rank-6) containing rep_j, members from orbits >= j: {m9}")
    say("        => every such 8-set lies inside orbit 11 (price 7), and orbit 11 has none of size 9: m = 9 impossible")
    R["m9_max_clique_by_orbit"] = m9

    if args.full:
        O11 = sum(1 << a for a, k in enumerate(V) if orb[k] == 11)
        r = vpos[reps[11]]; cnt = [0]
        def ext(P, need_):
            if need_ == 0:
                cnt[0] += 1; return
            while P:
                if P.bit_count() < need_:
                    return
                b = P & -P; i = b.bit_length() - 1; P ^= b
                ext(P & adj[i], need_ - 1)
        ext(adj[r] & O11, 7)
        say(f"--full: 8-sets of pairwise miss-disjoint orbit-11 supports containing rep_11: {cnt[0]:,}")
        R["orbit11_8sets_containing_rep"] = cnt[0]

    # ---- step 3: m = 10 -----------------------------------------------------------------------
    sets7 = []; by_orbit = {}
    def dfs(chosen, c1, c2, pool, needn, j):
        if needn == 0:
            sets7.append((tuple(chosen), c1, c2)); by_orbit[j] = by_orbit.get(j, 0) + 1; return
        for i, k in enumerate(pool):
            if len(pool) - i < needn:
                return
            mk = miss[k]
            if mk & c2:
                continue
            dfs(chosen + [k], c1 | mk, c2 | (c1 & mk), pool[i + 1:], needn - 1, j)
    for j, rk in enumerate(reps):
        if pr[rk] != 6:
            continue
        dfs([rk], miss[rk], 0, [k for k in R6 if orb[k] >= j and k != rk], 6, j)
    need(len(sets7) == 1248 and by_orbit == {7: 1248}, f"triple-free 7-sets {by_orbit}")
    say(f"step 3a (m=10): triple-free 7-sets of rank-6 supports up to the group: {len(sets7):,} (all with smallest orbit 7)")
    leaves = 0; covering = 0
    for chosen, c1, c2 in sets7:
        cov7 = 0
        for k in chosen: cov7 |= sup[k]
        pool = [k for k in range(K) if not miss[k] & c2]
        byp = {p: [k for k in pool if pr[k] == p] for p in (6, 7, 8, 9)}
        for p1, p2, p3 in TRIPLES:                       # price multisets with p1 <= p2 <= p3, sum <= 21
            A1 = byp[p1]
            for a, k1 in enumerate(A1):
                m1 = miss[k1]; d1 = c2 | (c1 & m1); e1 = c1 | m1; s1 = cov7 | sup[k1]
                A2 = [k for k in (A1[a:] if p2 == p1 else byp[p2]) if not miss[k] & d1]
                for b, k2 in enumerate(A2):
                    m2 = miss[k2]; d2 = d1 | (e1 & m2); s2 = s1 | sup[k2]
                    for k3 in (A2[b:] if p3 == p2 else byp[p3]):
                        if miss[k3] & d2:
                            continue
                        leaves += 1
                        if (s2 | sup[k3]) == LIVE:
                            covering += 1
    need(leaves == 612304, f"frame-feasible completions {leaves}")
    need(covering == 0, "a 10-term completion covers every line")
    say(f"step 3b: completions by 3 arbitrary supports (multiset, total price <= 21) with no frame missed 3x: "
        f"{leaves:,}; covering every line: {covering} => m = 10 impossible")
    R.update(sets7=len(sets7), m10_frame_feasible_completions=leaves, m10_covering=covering)
    say("THEOREM: no dictionary decoder of cost <= 63 (frames >= 8, lines >= 1, price >= 6); minimum 64")

    # ---- step 4: extension by atom43 ------------------------------------------------------------
    R.update(extension(G, sup, pr, FR, miss, hit, orb, reps, xc=False))
    return R


def load_atom(G):
    A = json.loads((HERE / "data" / "atom43.json").read_text())
    pos = {r: i for i, r in enumerate(G["roots"])}
    rows = [pos[tuple(x["root"])] for x in A["rows"]]; cols = [pos[tuple(x["root"])] for x in A["cols"]]
    need(all(i in G["live"] for i in rows + cols), "atom lines are live")
    S = [[1 if c["level"] <= r["level"] else 0 for c in A["cols"]] for r in A["rows"]]
    return rows, cols, S


def extension(G, sup, pr, FR, miss, hit, orb, reps, xc):
    R = {}
    rows, cols, S = load_atom(G)
    B = bd_block(G, rows, cols)
    MS = [[B[i][j] * S[i][j] for j in range(len(cols))] for i in range(len(rows))]
    cost = f2rank_low([sum(v << j for j, v in enumerate(r)) for r in MS])
    rq = qrank(MS)
    need(cost == 3 and rq == 4, f"atom43 cost {cost} rankQ {rq}")
    ra = []
    for f in FR:
        ri = [i for i, r in enumerate(rows) if f >> r & 1]; ci = [j for j, c in enumerate(cols) if f >> c & 1]
        sub = [[S[i][j] for j in ci] for i in ri]
        ra.append((f2rank_low([sum(v << t for t, v in enumerate(r)) for r in sub]) if xc else qrank(sub)) if ri and ci else 0)
    nA = sum(1 for x in ra if x >= 1)
    need(nA == 101 and max(ra) == 1, f"atom43 frame ranks: {nA} frames, max {max(ra)}")
    A = sum(1 << j for j, x in enumerate(ra) if x)
    frank = sum(ra)
    say(f"step 4 (atom43): cost 3, rankQ 4, rank 1 on each of {nA} frames, 0 elsewhere (frame-rank {frank})")
    K = len(sup); NF = len(FR)
    eff = max(Fraction(h.bit_count(), p) for h, p in zip(hit, pr))
    need(eff == 170, f"max dictionary frame efficiency {eff}")
    need(170 * (BUDGET - 3) + frank >= 8 * NF > 170 * (BUDGET - 6) + 2 * frank, "efficiency cap arithmetic")
    say("  4a: dictionary frames-met/price <= 170; 8*1260 = 10080 > 170*(63-3k) + 101k for k >= 2 => at most one atom")
    say("  4b: m_D <= 7 dead (a frame the atom misses needs 8 dictionary hits); m_D >= 11 costs >= 66")
    mm7 = min(miss[k].bit_count() for k in range(K) if pr[k] <= 7)
    need(2 * mm7 > nA and 6 + 7 * 8 > BUDGET - 3, "m_D = 8 counting")
    say(f"  4c: m_D = 8: dictionary misses only on atom frames, once each (<= {nA}); price<=7 supports miss >= {mm7} "
        f"frames => at most one => c_D >= 62 > 60")
    R6 = [k for k in range(K) if pr[k] == 6]
    notA = ((1 << NF) - 1) ^ A
    if xc:
        bad = sum(1 for a in range(len(R6)) for b in range(a, len(R6)) if not (miss[R6[a]] & miss[R6[b]] & notA))
    else:                                   # bitset form: b is bad for a iff b meets every frame in miss_a - A
        H6 = [0] * NF
        for a, k in enumerate(R6):
            for j in bits(hit[k]):
                H6[j] |= 1 << a
        bad = 0
        for a, k in enumerate(R6):
            x = (1 << len(R6)) - 1
            for j in bits(miss[k] & notA):
                x &= H6[j]
                if not x:
                    break
            bad += x.bit_count()
    need(bad == 0, "a rank-6 pair whose shared misses lie inside the atom frames")
    say(f"  4d: m_D = 9 (c_D <= 60 => >= 3 rank-6 terms): no pair of rank-6 supports (repeats included, "
        f"{len(R6) * (len(R6) + 1) // 2:,} pairs) has its shared missed frames inside the atom frames => dead")
    rq6 = max(qrank(bd_block(G, list(bits(sup[reps[o]])), list(bits(sup[reps[o]])))) for o in range(len(reps))
              if pr[reps[o]] == 6)
    need(rq6 == 7 and 10 * rq6 + rq < 112, "global rank row")
    say(f"  4e: m_D = 10 => ten rank-6 terms; rank over Q: 10*{rq6} + {rq} = {10 * rq6 + rq} < 112 => dead")
    say("EXTENSION: dictionary + atom43 orbit (either orientation) has no decoder of cost <= 63")
    R.update(atom_frames=nA, atom_frame_rank=frank, max_efficiency=170, min_miss_price_le7=mm7, rank6_rankQ=rq6)
    return R


# ----------------------------------------------------------------------------------------------
# independent cross-check (auditor's construction: row-span prices, direct frame scan, union-find
# orbits in a different labelling, colouring-bound clique search, line-driven completion search)
# ----------------------------------------------------------------------------------------------
def cross_check(args):
    R = {}
    mod, res5088 = load_pr5088(args.pr_dir, args.repo)
    roots, dots, live, pairing, adj = mod["geometry"]()
    G = own_geometry()
    # dictionary straight from #5088's code; prices re-derived by explicit row-span enumeration
    atoms = mod["dictionary"](dots, live, pairing)
    S = [s for s, _ in atoms]; P = []
    for s, r in atoms:
        span = {0}
        for i in bits(s):
            row = pairing[i] & s
            span |= {v ^ row for v in span}
        p = len(span).bit_length() - 1
        need(len(span) == 1 << p and p == r, "row-span price == #5088 rank")
        P.append(p)
    need(dict_hash(S, P) == DICT_SHA, "dictionary sha256")
    say("xc dictionary: #5088 dictionary(), every price re-derived by row-span enumeration; sha256 matches")
    L = list(bits(live))
    frames = []
    def rec(chosen, start):
        if len(chosen) == 8:
            frames.append(sum(1 << c for c in chosen)); return
        for t in range(start, len(L)):
            c = L[t]
            if all(dots[c][d] == 0 for d in chosen):
                rec(chosen + [c], t + 1)
    rec([], 0)
    need(len(frames) == 1260, "xc frames")
    frames.sort()
    miss = []
    for s in S:
        m = 0
        for k, f in enumerate(frames):
            if not s & f:
                m |= 1 << k
        miss.append(m)
    hit = [((1 << 1260) - 1) ^ m for m in miss]
    say("xc frames: 1,260 by pairwise-orthogonality recursion; miss sets by direct scan")
    gens = json.loads((HERE / "data" / "aut_generators.json").read_text())["generators"]
    for g in gens:
        need(all((dots[a][b] == 0) == (dots[g[a]][g[b]] == 0) for a in L for b in L), "xc automorphism")
    idx = {s: i for i, s in enumerate(S)}
    parent = list(range(len(S)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for g in gens:
        for i, s in enumerate(S):
            t = idx[permute(s, g)]
            need(P[t] == P[i], "price invariant")
            parent[find(i)] = find(t)
    groups = {}
    for i in range(len(S)):
        groups.setdefault(find(i), []).append(i)
    orbs = sorted(groups.values(), key=lambda o: (P[o[0]], len(o), S[o[0]].bit_count(), o[0]))
    oid = [0] * len(S)
    for k, o in enumerate(orbs):
        for i in o:
            oid[i] = k
    need(len(orbs) == 16, "xc 16 orbits")
    say("xc orbits (union-find, labels sorted by price/size): " + " ".join(f"{k}:{len(o)}/p{P[o[0]]}" for k, o in enumerate(orbs)))
    K = len(S)
    # m <= 8
    need(min(P) == 6 and all(P[i] >= 8 for i in range(K) if miss[i] == 0), "xc m<=8")
    # rank-6 pairwise (direct double loop)
    R6 = [i for i in range(K) if P[i] == 6]
    need(all(miss[R6[a]] & miss[R6[b]] for a in range(len(R6)) for b in range(a + 1, len(R6))), "xc rank-6 pairs")
    say(f"xc m<=8 ok; rank-6 pairs share a missed frame: all {len(R6) * (len(R6) - 1) // 2:,} (direct loop)")
    # m = 9: colouring-bound maximum clique through each rank-6/7 orbit representative
    V = [i for i in range(K) if P[i] in (6, 7)]; n = len(V)
    adjv = []
    for a in range(n):
        ma = miss[V[a]]
        adjv.append(sum(1 << b for b in range(n) if b != a and not ma & miss[V[b]]))
    six = [P[V[a]] == 6 for a in range(n)]
    def colour_bound(Pm):
        cnt = 0; U = Pm
        while U:
            cnt += 1; Q = U
            while Q:
                b = Q & -Q; i = b.bit_length() - 1; Q &= ~b & ~adjv[i]; U &= ~b
        return cnt
    def maxcl(r, Pm, cap):
        best = [1]
        def go(size, Pm, n6):
            best[0] = max(best[0], size)
            if best[0] >= cap or not Pm or size + colour_bound(Pm) <= best[0]:
                return
            while Pm:
                if size + Pm.bit_count() <= best[0]:
                    return
                b = Pm & -Pm; i = b.bit_length() - 1; Pm ^= b
                if n6 + six[i] > 1:
                    continue
                go(size + 1, Pm & adjv[i], n6 + six[i])
        go(1, Pm, six[r])
        return best[0]
    m9 = {}
    pos = {k: a for a, k in enumerate(V)}
    for j, o in enumerate(orbs):
        if P[o[0]] not in (6, 7):
            continue
        r = pos[o[0]]
        Pm = adjv[r] & sum(1 << a for a in range(n) if oid[V[a]] >= j)
        m9[j] = maxcl(r, Pm, 9)
    # the only orbit allowed to reach 8 is the price-7 orbit of size 448 with 28-line supports
    big = [j for j, v in m9.items() if v >= 8]
    need(all(v <= 8 for v in m9.values()) and len(big) == 1 and P[orbs[big[0]][0]] == 7 and len(orbs[big[0]]) == 448
         and S[orbs[big[0]][0]].bit_count() == 28, f"xc m=9 {m9}")
    # every 8-set containing rep of that orbit stays inside it: m9[j] for its orbit computed among orbits >= j,
    # so also require: the max set containing rep_j for all other j < 8
    say(f"xc m=9: max pairwise miss-disjoint set (colouring bound) by orbit label {m9}: only the 28-line price-7 orbit reaches 8")
    # m = 10: triple-free 7-sets in this labelling + line-driven completion (budget 63 - 42 = 21)
    sets7 = []                     # bitset DFS: candidates that would miss a doubly-missed frame are dropped
    n6 = len(R6); HIT6 = [0] * 1260
    for a, i in enumerate(R6):
        for f in range(1260):
            if not miss[i] >> f & 1:
                HIT6[f] |= 1 << a
    for j, o in enumerate(orbs):
        if P[o[0]] != 6:
            continue
        rep = o[0]; ra = R6.index(rep)
        C0 = sum(1 << a for a, i in enumerate(R6) if oid[i] >= j and i != rep)
        def dfs(cl, M1, M2, C):
            if len(cl) == 7:
                sets7.append(([R6[a] for a in cl], M1, M2)); return
            while C:
                if len(cl) + C.bit_count() < 7:
                    return
                b = C & -C; a = b.bit_length() - 1; C ^= b
                m = miss[R6[a]]; new2 = M1 & m & ~M2; C2 = C
                for f in bits(new2):
                    C2 &= HIT6[f]
                    if not C2:
                        break
                dfs(cl + [a], M1 | m, M2 | new2, C2)
        dfs([ra], miss[rep], 0, C0)
    need(len(sets7) == 2208, f"xc triple-free 7-sets {len(sets7)} (auditor's count: 1,536 + 672)")
    say(f"xc m=10: {len(sets7):,} triple-free 7-sets of rank-6 supports in this labelling (= auditor's 1,536 + 672)")
    order = sorted(range(K), key=lambda i: P[i])
    def complete(left, bud, M1, M2, U, okl):
        """line-driven: some remaining term covers the lowest uncovered line."""
        if left == 0:
            return U == 0
        if bud < 6 * left:
            return False
        u = (U & -U).bit_length() - 1 if U else None
        for t in okl:
            if P[t] > bud - 6 * (left - 1):
                break
            if u is not None and not S[t] >> u & 1:
                continue
            if left == 1 and (S[t] & U) != U:
                continue
            m = miss[t]
            if m & M2:
                continue
            if complete(left - 1, bud - P[t], M1 | m, M2 | (M1 & m), U & ~S[t], okl):
                return True
        return False
    def start(cl, M1, M2, bud):
        cov = 0
        for c in cl: cov |= S[c]
        okl = [t for t in order if not miss[t] & M2]
        return complete(3, bud, M1, M2, live & ~cov, okl)
    sols = sum(1 for cl, M1, M2 in sets7 if start(cl, M1, M2, BUDGET - 42))
    need(sols == 0, "xc m=10 completion found")
    ctl = next((1 for cl, M1, M2 in sets7[:200] if start(cl, M1, M2, 45)), 0)
    need(ctl == 1, "xc non-vacuity control: a 45-unit completion budget finds a cover")
    say(f"xc m=10: line-driven completion at budget 21 covers no line set: 0 of {len(sets7):,} "
        f"(control: budget 45 finds a cover)")
    R.update(xc_sets7=len(sets7), xc_m9=m9)
    need(roots == G["roots"], "xc geometry")
    H = dict(G)
    R.update(extension(H, S, P, frames, miss, hit, None, [o[0] for o in orbs], xc=True))
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pr-dir", help="directory holding #5088's verify.py and RESULT.json (default: vendor/ copy)")
    ap.add_argument("--repo", help="read #5088's files with git show from this checkout (needs the PR head fetched)")
    ap.add_argument("--cross-check", action="store_true", help="run the independent second implementation")
    ap.add_argument("--full", action="store_true", help="also count orbit-11 8-sets (about 2.5 minutes)")
    ap.add_argument("--output", help="write a JSON receipt to this new path")
    args = ap.parse_args()
    mode = "cross-check" if args.cross_check else "primary"
    try:
        R = cross_check(args) if args.cross_check else primary(args)
        status = "PASS"
    except Exception as e:  # noqa: BLE001
        say(f"FAILED: {type(e).__name__}: {e}")
        R, status = {}, "FAIL"
    R.update(mode=mode, status=status, pr5088_head=PR_HEAD, dictionary_sha256=DICT_SHA,
             seconds=round(time.time() - T0, 1), python=sys.version.split()[0])
    if args.output:
        with open(args.output, "x") as fh:
            json.dump(R, fh, indent=1, sort_keys=True, default=str); fh.write("\n")
    print(f"{status}_P112_DICTIONARY_FLOOR_64 ({mode}, {R['seconds']} s)")
    sys.exit(0 if status == "PASS" else 1)


if __name__ == "__main__":
    main()
