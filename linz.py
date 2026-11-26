"""Linearized construction for ISA-L points (A = {0..n-1} inside the GF(2)-subspace V = {0..15}), r = 4.
Check polynomials g(y) = a*y^3 + b*y + c  (y = x - a*):  g(y)*y = a*y^4 + b*y^2 + c*y is GF(2)-linear in y.
For a 4-dim GF(2)-subspace Y of GF(256):  S_Y = { (a,b,c) : a*d^4 + b*d^2 + c*d in Y for all d in V }.
Every g in S_Y has g(d) in Y/d, so each helper sends at most 4 bits."""
import sys, itertools, json
sys.path.insert(0, '/home/claude/rs')
from verify import mul, power, inv, rank_gf2, verify

def subspaces_dim4():
    """all 4-dim subspaces of GF(2)^8 as reduced row echelon bases (pivot columns distinct, descending)."""
    out = []
    for pivs in itertools.combinations(range(7, -1, -1), 4):          # pivot bit positions, descending
        free = [b for b in range(8) if b not in pivs]
        # row t has pivot pivs[t]; free bits below its pivot are arbitrary
        choices = []
        for p in pivs:
            fb = [b for b in free if b < p]
            choices.append([(1 << p) | sum((1 << fb[i]) for i in range(len(fb)) if m >> i & 1) for m in range(1 << len(fb))])
        for rows in itertools.product(*choices):
            out.append(rows)
    return out

def annihilator(rows):
    """4 functionals (bitmasks f; value = parity(f & v)) vanishing on span(rows)."""
    # solve f . r = 0 for all rows: nullspace of 4x8 matrix
    sols = []
    for f in range(1, 256):
        if all(bin(f & r).count('1') % 2 == 0 for r in rows):
            sols.append(f)
    basis = []
    for f in sols:
        if rank_gf2(basis + [f]) > len(basis):
            basis.append(f)
        if len(basis) == 4:
            break
    return basis

def nullspace(masks, nbits):
    """nullspace (list of bit-vectors) of GF(2) system: each mask is a row over nbits unknowns."""
    rows = []; pivcols = []
    for m in masks:
        for r, pc in zip(rows, pivcols):
            if m >> pc & 1: m ^= r
        if m:
            pc = m.bit_length() - 1
            for i in range(len(rows)):
                if rows[i] >> pc & 1: rows[i] ^= m
            rows.append(m); pivcols.append(pc)
    free = [c for c in range(nbits) if c not in pivcols]
    ns = []
    for fcol in free:
        v = 1 << fcol
        for r, pc in zip(rows, pivcols):
            if r >> fcol & 1: v |= 1 << pc
        ns.append(v)
    return ns

Vbasis = [1, 2, 4, 8]
# contribution of unknown bit (coef index t in {0:a,1:b,2:c}, bit position q) to l(d) for d in Vbasis
contrib = {}
for d in Vbasis:
    pw = [power(d, 4), power(d, 2), d]
    for t in range(3):
        for q in range(8):
            contrib[(d, t * 8 + q)] = mul(1 << q, pw[t])

def S_of(Yrows):
    ann = annihilator(Yrows)
    masks = []
    for d in Vbasis:
        for f in ann:
            m = 0
            for bit in range(24):
                if bin(f & contrib[(d, bit)]).count('1') & 1: m |= 1 << bit
            masks.append(m)
    return nullspace(masks, 24)

def abc(v): return (v & 255, (v >> 8) & 255, (v >> 16) & 255)     # a, b, c

def ranks_for(S, deltas):
    out = []
    for d in deltas:
        d4, d2 = power(d, 4), power(d, 2); di = inv(d)
        vals = []
        for v in S:
            a, b, c = abc(v)
            vals.append(mul(di, mul(a, d4) ^ mul(b, d2) ^ mul(c, d)))
        out.append(rank_gf2(vals))
    return out

if __name__ == '__main__':
    subs = subspaces_dim4()
    print(len(subs), "four-dimensional subspaces Y")
    stats = {}; best = None
    for Y in subs:
        S = S_of(Y)
        if len(S) != 8:
            stats[('dim', len(S))] = stats.get(('dim', len(S)), 0) + 1
            continue
        if rank_gf2([abc(v)[2] for v in S]) != 8:
            stats['c-not-injective'] = stats.get('c-not-injective', 0) + 1
            continue
        rk = ranks_for(S, range(1, 16))
        key = tuple(sorted(rk)); stats[key] = stats.get(key, 0) + 1
        tot = sum(rk)
        if best is None or tot < best[0]:
            best = (tot, Y, S, rk)
    for k, v in sorted(stats.items(), key=lambda kv: str(kv[0])): print(k, v)
    print("best total for n=16:", best[0], "ranks", best[3], "Y basis", best[1])
    json.dump({"Y": list(best[1]), "S": [list(abc(v)) for v in best[2]]}, open('/tmp/linz_best.json', 'w'))
