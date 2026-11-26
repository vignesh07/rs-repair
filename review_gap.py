#!/usr/bin/env python3
"""Independent review of the RS(6,2) point-set gap (RS6_2_POINTSET_GAP.md), from the check-polynomial side.

It shares only GF(256) arithmetic and the repair verifier (verify.py) with the note under review; the argument is different.

Dual-side derivation. A linear scheme for failed point t is an 8-dimensional GF(2)-space D of polynomials of degree < 4
with evaluation at t injective; helper j sends dim{g(j) : g in D} bits. Evaluation at any four helpers is injective on
D (a cubic with four roots is zero), so the bits of any four helpers add up to at least 8, and a 10-bit scheme has
exactly 2 bits at each of the five helpers. Then evaluation at helpers 1..4 identifies D with V1 x V2 x V3 x V4
(2-dimensional value spaces), i.e. D = sum_i V_i * l_i(x) with l_i the Lagrange basis on those four points. Two bits at
the fifth helper force l_i(a5) V_i to be one and the same space W for all i. So every 10-bit scheme is
    D(W) = sum_i (W / l_i(a5)) * l_i(x)     for a 2-dimensional W,
and it repairs t exactly when sum_i mu_i W = GF(256), mu_i = l_i(t) / l_i(a5).
Here mu_i is a common constant times the note's c_i = (i - h)/(i - t) with h = a5.

This script (1) enumerates all 10,795 spaces W for every node of both point sets and every choice of fifth helper,
(2) hands one D(W) per repairable node to the real-repair verifier, (3) recounts the six-subsets of GF(16)."""
import itertools, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify import mul, inv, power, rank_gf2, verify


def polymul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] ^= mul(a, b)
    return out


def lagrange(pts4):
    """coefficient lists (lowest power first) of l_i, i = 0..3"""
    out = []
    for i, a in enumerate(pts4):
        num, den = [1], 1
        for j, b in enumerate(pts4):
            if j != i:
                num = polymul(num, [b, 1]); den = mul(den, a ^ b)
        out.append([mul(c, inv(den)) for c in num])
    return out


def ev(p, x):
    acc = 0
    for c in reversed(p):
        acc = mul(acc, x) ^ c
    return acc


def two_dim_subspaces():
    seen = set()
    for a in range(1, 256):
        for b in range(a + 1, 256):
            key = tuple(sorted((a, b, a ^ b)))
            if key not in seen:
                seen.add(key); yield (a, b)


SUBS = list(two_dim_subspaces())
assert len(SUBS) == 10795


def analyse(points, t, a5):
    """rank of the mu_i, number of W that give a 10-bit repair, and one such scheme"""
    others = [p for p in points if p not in (t, a5)]
    L = lagrange(others)
    lam = [ev(l, a5) for l in L]; mu = [mul(ev(l, t), inv(x)) for l, x in zip(L, lam)]
    good, example = 0, None
    for w1, w2 in SUBS:
        if rank_gf2([mul(m, w) for m in mu for w in (w1, w2)]) == 8:
            good += 1
            if example is None:
                example = [[mul(mul(w, inv(lam[i])), c) for c in L[i]] for i in range(4) for w in (w1, w2)]
    return rank_gf2(mu), good, example


def profile(label, points, code):
    print(label, points)
    for t in points:
        ranks, goods, example = [], [], None
        for a5 in points:
            if a5 == t:
                continue
            r, g, ex = analyse(points, t, a5)
            ranks.append(r); goods.append(g); example = example or ex
        assert len(set(ranks)) == 1 and len(set(goods)) == 1, "the answer must not depend on which helper is called the fifth"
        line = f"  node at point {t:3d}: rank of the four multipliers {ranks[0]}, spaces W giving a 10-bit repair: {goods[0]} of 10795"
        if example:
            bw, rk = verify({"n": 6, "k": 2, "points": list(points), "node": list(points).index(t), "polys": example, "code": code}, trials=2000, quiet=True)
            assert bw == 10
            line += f"; one verified with 2000 real repairs, helper bits {[rk[j] for j in sorted(rk)]}"
        else:
            line += "; no 10-bit linear repair exists"
        print(line)


if __name__ == "__main__":
    profile("ISA-L points", (0, 1, 2, 3, 4, 5), "isal")
    w = power(2, 17)
    profile("subfield points as printed in the thesis appendix", (0, 1, w, power(w, 2), power(w, 3), power(w, 4)), "eval")
    z = 78
    profile("subfield points as written in the note under review", (0, 1, z, power(z, 2), power(z, 3), power(z, 4)), "eval")
    f16 = [x for x in range(256) if power(x, 16) == x]
    hist = {}
    for pts in itertools.combinations(f16, 6):
        bad = 0
        for t in pts:
            a5 = next(p for p in pts if p != t)
            others = [p for p in pts if p not in (t, a5)]
            # mu_i up to a common constant: (a_i - a5)/(a_i - t)
            if rank_gf2([mul(a ^ a5, inv(a ^ t)) for a in others]) < 4:
                bad += 1
        hist[bad] = hist.get(bad, 0) + 1
    print("six-subsets of GF(16) by number of nodes that cannot be repaired with 10 bits:", dict(sorted(hist.items())))
