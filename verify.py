#!/usr/bin/env python3
"""Independent end-to-end verifier for linear (trace) repair schemes of Reed-Solomon codes over GF(256).

It shares no code with the search programs. For a claimed scheme it
  1. builds the code's generator matrix (ISA-L Cauchy construction, or a generic GRS code on given points),
  2. checks that the 8 check polynomials give dual codewords,
  3. recomputes the bandwidth (sum of GF(2)-ranks at the helpers),
  4. actually repairs the lost symbol of random codewords using only those bits.

Scheme format (JSON): {"n":..,"k":..,"points":[..],"node":s,"polys":[[c0,c1,..,c_{r-1}], ... 8 polys]}
  polys are coefficient lists in powers of x (lowest first), entries are bytes in GF(256) with polynomial 0x11D.
"""
import json, random, sys

POLY = 0x11D


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        if a & 0x100:
            a ^= POLY
        b >>= 1
    return r


def power(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def inv(a):
    assert a != 0
    return power(a, 254)


def trace(a):
    t, x = 0, a
    for _ in range(8):
        t ^= x
        x = mul(x, x)
    assert t in (0, 1)
    return t


def rank_gf2(vals):
    basis = []
    for v in vals:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def basis_and_coords(vals):
    """greedy GF(2)-basis of the span of vals, and for each val its expansion as a subset (bitmask) of the basis."""
    basis, red, combo = [], [], []          # red[i] = reduced vector, combo[i] = mask over basis indices
    coords = []
    for v in vals:
        x, m = v, 0
        for i, rv in enumerate(red):
            if x ^ rv < x:
                x ^= rv; m ^= combo[i]
        if x:
            basis.append(v); red.append(x); combo.append(m ^ (1 << (len(basis) - 1)))
            coords.append(1 << (len(basis) - 1))
        else:
            coords.append(m)
    # sanity: every val equals XOR of basis elements selected by coords
    for v, m in zip(vals, coords):
        acc = 0
        for i, b in enumerate(basis):
            if m >> i & 1:
                acc ^= b
        assert acc == v
    return basis, coords


def poly_eval(coeffs, x):
    acc = 0
    for c in reversed(coeffs):
        acc = mul(acc, x) ^ c
    return acc


def isal_generator(n, k):
    """ISA-L gf_gen_cauchy1_matrix: identity on top, then a[i][j] = inv(i ^ j) for i >= k."""
    G = [[1 if i == j else 0 for j in range(k)] for i in range(k)]
    for i in range(k, n):
        G.append([inv(i ^ j) for j in range(k)])
    return G


def gf_matinv(M):
    n = len(M)
    A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        p = next(i for i in range(col, n) if A[i][col])
        A[col], A[p] = A[p], A[col]
        iv = inv(A[col][col])
        A[col] = [mul(iv, x) for x in A[col]]
        for i in range(n):
            if i != col and A[i][col]:
                f = A[i][col]
                A[i] = [x ^ mul(f, y) for x, y in zip(A[i], A[col])]
    return [row[n:] for row in A]


def vandermonde_generator(n, k):
    """Backblaze JavaReedSolomon / klauspost-reedsolomon default: Vandermonde rows i^j (i = 0..n-1), made systematic."""
    V = [[power(i, j) if j else 1 for j in range(k)] for i in range(n)]
    top_inv = gf_matinv(V[:k])
    G = []
    for i in range(n):
        G.append([0] * k)
        for j in range(k):
            acc = 0
            for t in range(k):
                acc ^= mul(V[i][t], top_inv[t][j])
            G[i][j] = acc
    assert all(G[i][j] == (1 if i == j else 0) for i in range(k) for j in range(k))
    return G


def grs_dual_multipliers(G, points, r):
    """Find column multipliers lam with sum_node lam[node]*g(points[node])*c[node] = 0 for all codewords c, deg g < r.
    Solved numerically: lam is a nonzero vector with (lam_i * a_i^m)_i in the dual for m = 0..r-1."""
    n, k = len(G), len(G[0])
    # unknown lam (n bytes); constraints: for each data column d and each m: sum_i lam_i a_i^m G[i][d] = 0  (GF(256)-linear)
    rows = []
    for m in range(r):
        for d in range(k):
            rows.append([mul(power(points[i], m) if m else 1, G[i][d]) for i in range(n)])
    # Gaussian elimination over GF(256) to get the null space
    M = [row[:] for row in rows]
    piv = []
    rr = 0
    for col in range(n):
        p = next((i for i in range(rr, len(M)) if M[i][col]), None)
        if p is None:
            continue
        M[rr], M[p] = M[p], M[rr]
        iv = inv(M[rr][col])
        M[rr] = [mul(iv, x) for x in M[rr]]
        for i in range(len(M)):
            if i != rr and M[i][col]:
                f = M[i][col]
                M[i] = [x ^ mul(f, y) for x, y in zip(M[i], M[rr])]
        piv.append(col)
        rr += 1
    free = [c for c in range(n) if c not in piv]
    assert len(free) == 1, f"expected a 1-dim solution space for the multipliers, got {len(free)}"
    lam = [0] * n
    lam[free[0]] = 1
    for i, pc in enumerate(piv):
        lam[pc] = M[i][free[0]]
    assert all(lam), "multipliers must be nonzero for an MDS/GRS code"
    return lam


def verify(scheme, trials=2000, seed=1, quiet=False):
    n, k, pts, s = scheme["n"], scheme["k"], scheme["points"], scheme["node"]
    r = n - k
    polys = scheme["polys"]
    assert len(polys) == 8 and all(len(p) <= r for p in polys), "need 8 polynomials of degree < r"
    assert len(set(pts)) == n
    code = scheme.get("code", "isal")
    if code == "isal":
        assert pts == list(range(n)), "ISA-L code uses points 0..n-1"
        G = isal_generator(n, k)
    elif code == "vandermonde":
        assert pts == list(range(n)), "Backblaze-style code uses points 0..n-1"
        G = vandermonde_generator(n, k)
    elif code == "eval":
        # plain evaluation code on arbitrary distinct points: codeword = (f(p) for p in points), deg f < k
        G = [[power(p, j) if j else 1 for j in range(k)] for p in pts]
    else:
        raise ValueError("unknown code")
    lam = grs_dual_multipliers(G, pts, r)
    # 1) each check polynomial gives a dual codeword (direct orthogonality test against the generator matrix)
    for g in polys:
        h = [mul(lam[i], poly_eval(g, pts[i])) for i in range(n)]
        for d in range(k):
            acc = 0
            for i in range(n):
                acc ^= mul(h[i], G[i][d])
            assert acc == 0, "check polynomial is not in the dual code"
    # 2) full rank at the failed node, ranks at helpers
    theta = [mul(lam[s], poly_eval(g, pts[s])) for g in polys]
    assert rank_gf2(theta) == 8, "values at the failed node do not span GF(256)"
    helpers = [j for j in range(n) if j != s]
    gam = {j: [mul(lam[j], poly_eval(g, pts[j])) for g in polys] for j in helpers}
    ranks = {j: rank_gf2(gam[j]) for j in helpers}
    bw = sum(ranks.values())
    bases = {j: basis_and_coords(gam[j]) for j in helpers}
    assert all(len(bases[j][0]) == ranks[j] for j in helpers)
    # 3) real repairs
    rng = random.Random(seed)
    for _ in range(trials):
        data = [rng.randrange(256) for _ in range(k)]
        cw = []
        for i in range(n):
            acc = 0
            for d in range(k):
                acc ^= mul(G[i][d], data[d])
            cw.append(acc)
        sent = 0
        t = [0] * 8
        for j in helpers:
            basis, coords = bases[j]
            bits = [trace(mul(b, cw[j])) for b in basis]      # what helper j transmits
            sent += len(bits)
            for i in range(8):
                v = 0
                for q, bit in enumerate(bits):
                    if coords[i] >> q & 1:
                        v ^= bit
                t[i] ^= v
        assert sent == bw
        cands = [x for x in range(256) if all(trace(mul(theta[i], x)) == t[i] for i in range(8))]
        assert cands == [cw[s]], "repair failed"
    if not quiet:
        where = f"points 0..{n-1}" if pts == list(range(n)) else f"points {pts}"
        print(f"OK  RS({n},{k}) {code} code, {where}, node {s}: bandwidth {bw} bits (classical {8*k}), "
              f"helper ranks {[ranks[j] for j in helpers]}, {trials} random repairs correct")
    return bw, ranks


def scheme_from_c(n, r, node, cflat, points=None):
    """Convert the search program's output (c[i][m] = coefficient of (x-a*)^m in g_{e_i}, i=0..7, m=1..r-1) into plain polys."""
    pts = points or list(range(n))
    a = pts[node]
    polys = []
    it = iter(cflat)
    for i in range(8):
        cm = [1 << i] + [next(it) for _ in range(r - 1)]       # coefficients in powers of (x + a)
        # expand sum_m cm[m] (x+a)^m
        coeffs = [0] * r
        binom_pow = [1]                                           # (x+a)^0
        for m in range(r):
            for d, bc in enumerate(binom_pow):
                coeffs[d] ^= mul(cm[m], bc)
            # multiply binom_pow by (x + a)
            nxt = [0] * (len(binom_pow) + 1)
            for d, bc in enumerate(binom_pow):
                nxt[d + 1] ^= bc
                nxt[d] ^= mul(bc, a)
            binom_pow = nxt
        polys.append(coeffs)
    return {"n": n, "k": n - r, "points": pts, "node": node, "polys": polys, "code": "isal"}


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "selftest":
        # classical repair written as a trace scheme must cost 8*k bits: g_i = e_i * prod_{j in S}(x - a_j)/(same at a*), |S| = r-1
        n, k, s = 9, 6, 0
        r = n - k
        pts = list(range(n))
        S = [1, 2]
        base = [1]
        for a in S:
            nxt = [0] * (len(base) + 1)
            for d, c in enumerate(base):
                nxt[d + 1] ^= c
                nxt[d] ^= mul(c, pts[a])
            base = nxt
        norm = inv(poly_eval(base, pts[s]))
        polys = [[mul(mul(1 << i, norm), c) for c in base] for i in range(8)]
        bw, _ = verify({"n": n, "k": k, "points": pts, "node": s, "polys": polys, "code": "isal"})
        assert bw == 8 * k
        print("selftest passed: classical repair reproduced as a 48-bit scheme")
    else:
        for path in sys.argv[1:]:
            verify(json.load(open(path)))
