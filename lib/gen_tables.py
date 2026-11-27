#!/usr/bin/env python3
"""Turn the best-known trace-repair schemes for the points 0..n-1 into the table the C library loads (rsr_tables.h).

Per code RS(n,k) and failed node s the table holds
    fin[8]                      the lost byte is the XOR of fin[b] over the set bits b of an intermediate byte z
    per helper i:  r            bits helper i sends per chunk byte
                   hmask[r]     bit j it sends = parity(hmask[j] & chunk_byte)
                   dmask[r]     that bit is XORed into the bits of z selected by dmask[j]
for the plain evaluation code c_i = f(i), deg f < k, GF(256) with polynomial 0x11D (Backblaze JavaReedSolomon and
klauspost/reedsolomon default).  Codes that differ by a nonzero factor per node (Intel ISA-L Cauchy matrix) are handled
in the library by rescaling hmask and fin.  The intermediate byte z and the helper bases are chosen (by ./opt) to make
the number of XORs in the rebuild as small as we could find.

usage: cc -O2 -o opt opt.c && python3 gen_tables.py [dir with the scheme json files, default ..] > rsr_tables.h
       (takes about half a minute and self-checks every scheme on both codes)
"""
import json, os, random, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..")
sys.path.insert(0, SRC)                     # verify.py lives next to the scheme files
from verify import mul, inv, trace, poly_eval, basis_and_coords, isal_generator, vandermonde_generator, rank_gf2
par = lambda x: bin(x).count("1") & 1


def load_best():
    best = {}
    for fn in ("schemes_isal.json", "schemes_rs6_2_optimal.json", "schemes_rs20_17.json"):
        for name, s in json.load(open(os.path.join(SRC, fn))).items():
            assert s["points"] == list(range(s["n"]))
            key = (s["n"], s["k"], s["node"])
            if key not in best or s["bandwidth_bits"] < best[key]["bandwidth_bits"]:
                best[key] = s
    return best


def dual_basis(theta):
    M = [sum(trace(mul(theta[t], 1 << e)) << e for e in range(8)) for t in range(8)]
    delta = []
    for u in range(8):
        sol = [x for x in range(256) if all(par(M[t] & x) == (t == u) for t in range(8))]
        assert len(sol) == 1
        delta.append(sol[0])
    return delta


def plain_dual_multipliers(n):
    lam = []
    for i in range(n):
        p = 1
        for j in range(n):
            if j != i:
                p = mul(p, i ^ j)
        lam.append(inv(p))
    return lam


def cauchy_multipliers(n, k):
    """u with: (ISA-L Cauchy codeword)[i] = u[i] * f(i) for some f of degree < k."""
    u = []
    for i in range(n):
        p = 1
        for t in range(k):
            if t != i:
                p = mul(p, i ^ t)
        u.append(inv(p))
    return u


def raw_tables(s):
    """per helper: (beta, D) with lost = XOR_j D[j] * Tr(beta[j] * c_i)."""
    n, k, node, polys = s["n"], s["k"], s["node"], s["polys"]
    lam = plain_dual_multipliers(n)
    theta = [mul(lam[node], poly_eval(g, node)) for g in polys]
    assert rank_gf2(theta) == 8
    delta = dual_basis(theta)
    out = {}
    for i in range(n):
        if i == node:
            continue
        gam = [mul(lam[i], poly_eval(g, i)) for g in polys]
        basis, coords = basis_and_coords(gam)
        dcol = []
        for j in range(len(basis)):
            d = 0
            for t in range(8):
                if coords[t] >> j & 1:
                    d ^= delta[t]
            dcol.append(d)
        assert rank_gf2(dcol) == len(dcol)
        out[i] = (basis, dcol)
    assert sum(len(b) for b, _ in out.values()) == s["bandwidth_bits"]
    return out


def optimise(raw, A):
    """A = 8 row masks (z_b = parity(A[b] & lost)).  Returns fin, {helper: (hmask, dmask)}, edges."""
    TA = [sum(par(A[b] & x) << b for b in range(8)) for x in range(256)]
    assert len(set(TA)) == 256
    fin = [TA.index(1 << b) for b in range(8)]
    res, edges = {}, 0
    for i, (beta, D) in raw.items():
        r = len(D)
        cands = []
        for m in range(1, 1 << r):
            v = 0
            for j in range(r):
                if m >> j & 1:
                    v ^= D[j]
            cands.append((bin(TA[v]).count("1"), TA[v], v))
        cands.sort()
        piv, newD, dmask = {}, [], []
        for w, tv, v in cands:
            x = tv
            for b in range(7, -1, -1):
                if x >> b & 1:
                    if b in piv:
                        x ^= piv[b]
                    else:
                        piv[b] = x; newD.append(v); dmask.append(tv); x = 0
                        break
            if len(newD) == r:
                break
        assert len(newD) == r
        # functional masks in the new basis: coordinates of F(c) = XOR_j D[j]*Tr(beta[j]*c) along newD
        combo = {}
        for m in range(1 << r):
            v = 0
            for j in range(r):
                if m >> j & 1:
                    v ^= newD[j]
            combo[v] = m
        assert len(combo) == 1 << r
        hmask = [0] * r
        for e in range(8):
            f = 0
            for b, d in zip(beta, D):
                if trace(mul(b, 1 << e)):
                    f ^= d
            m = combo[f]
            for j in range(r):
                if m >> j & 1:
                    hmask[j] |= 1 << e
        res[i] = (hmask, dmask)
        edges += sum(bin(x).count("1") for x in dmask)
    return fin, res, edges


def selfcheck(s, fin, tab, rng, trials=40):
    """repair random codewords of both real-library codes using only the table, the way the C library does."""
    n, k, node = s["n"], s["k"], s["node"]
    for G, u in ((vandermonde_generator(n, k), [1] * n), (isal_generator(n, k), cauchy_multipliers(n, k))):
        hm = {i: [sum(par(h & mul(inv(u[i]), 1 << e)) << e for e in range(8)) for h in hmask] for i, (hmask, _) in tab.items()}
        fn = [mul(u[node], f) for f in fin]
        for _ in range(trials):
            data = [rng.randrange(256) for _ in range(k)]
            cw = []
            for i in range(n):
                acc = 0
                for d in range(k):
                    acc ^= mul(G[i][d], data[d])
                cw.append(acc)
            z = 0
            for i, (_, dmask) in tab.items():
                for h, d in zip(hm[i], dmask):
                    if par(h & cw[i]):
                        z ^= d
            lost = 0
            for b in range(8):
                if z >> b & 1:
                    lost ^= fn[b]
            assert lost == cw[node], (n, k, node)


best = load_best()
codes = sorted({(n, k) for n, k, _ in best}, key=lambda c: (c[0] - c[1], c[0]))
keys = [(n, k, node) for n, k in codes for node in range(n)]
raws = {key: raw_tables(best[key]) for key in keys}
inp = [str(len(keys))]
for key in keys:
    inp.append(str(len(raws[key])))
    for i in sorted(raws[key]):
        inp.append(" ".join(map(str, [len(raws[key][i][1])] + raws[key][i][1])))
opt = subprocess.run([os.path.join(HERE, "opt")], input="\n".join(inp) + "\n", capture_output=True, text=True, check=True)
As = [list(map(int, line.split()))[1:] for line in opt.stdout.strip().split("\n")]
assert len(As) == len(keys)

rng = random.Random(7)
blob, node_off, index, stats = [], [], [], {}
for key, A in zip(keys, As):
    n, k, node = key
    if node == 0:
        index.append((n, k, len(node_off)))
    fin, tab, edges = optimise(raws[key], A)
    _, _, edges_plain = optimise(raws[key], [1 << b for b in range(8)])
    selfcheck(best[key], fin, tab, rng)
    stats.setdefault((n, k), []).append((best[key]["bandwidth_bits"], edges, edges_plain))
    node_off.append(len(blob))
    blob += fin
    for i in range(n):
        if i == node:
            blob.append(0xFF)
            continue
        hmask, dmask = tab[i]
        blob += [len(hmask)] + hmask + dmask

w = sys.stdout.write
w("/* Generated by gen_tables.py from the scheme files of github.com/vignesh07/rs-repair. Do not edit. */\n")
w(f"#define RSR_NCODES {len(index)}\n")
w("static const struct { uint8_t n, k; uint16_t first; } rsr_codes[RSR_NCODES] = {\n")
w("".join(f"    {{{n}, {k}, {first}}},\n" for n, k, first in index) + "};\n")
w(f"static const uint32_t rsr_node_off[{len(node_off)}] = {{\n")
for i in range(0, len(node_off), 12):
    w("    " + ", ".join(map(str, node_off[i:i + 12])) + ",\n")
w("};\n")
w(f"static const uint8_t rsr_blob[{len(blob)}] = {{\n")
for i in range(0, len(blob), 24):
    w("    " + ",".join(map(str, blob[i:i + 24])) + ",\n")
w("};\n")
e = sys.stderr.write
e(f"{len(index)} codes, {len(node_off)} node schemes, {len(blob)} table bytes; self-check passed on both codes\n")
e("code        bits/byte  payload planes read per rebuilt 512 bytes (avg over lost chunks): optimised / before\n")
for (n, k), v in stats.items():
    if (n, k) in ((14, 10), (9, 6), (6, 3), (12, 8), (16, 12), (8, 4), (6, 4), (10, 8), (5, 3), (4, 2), (20, 17)):
        e(f"RS({n},{k})".ljust(12) + f"{sum(x[0] for x in v)/len(v):7.1f} {sum(x[1] for x in v)/len(v):12.1f} / {sum(x[2] for x in v)/len(v):.1f}\n")
