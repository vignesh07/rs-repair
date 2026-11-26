"""r = 4, ISA-L points: best scheme per (n, node) inside the linearized family D = S_Y (and hyperplanes of the 9-dim S_Y)."""
import sys, json, itertools, pickle
sys.path.insert(0, '/home/claude/rs')
from verify import mul, power, inv, rank_gf2
from linz import subspaces_dim4, S_of, abc, ranks_for

def all_hyperplanes(S9):
    """8-dim subspaces of a 9-dim space given by basis S9 (as kernels of the 511 nonzero functionals)."""
    out = []
    for f in range(1, 512):
        # elements sum_{i in mask} S9[i] with parity(mask & f) == 0: basis = pick pivot p in f, then e_i (+ e_p if f_i)
        p = f.bit_length() - 1
        basis = []
        for i in range(9):
            if i == p: continue
            v = S9[i]
            if f >> i & 1: v ^= S9[p]
            basis.append(v)
        out.append(basis)
    return out

cands = {}       # rank vector (delta=1..15) -> D (list of 8 ints encoding (a,b,c))
subs = subspaces_dim4()
seen_special = 0
for idx, Y in enumerate(subs):
    S = S_of(Y)
    Ds = []
    if len(S) == 8: Ds = [S]
    elif len(S) == 9 and seen_special < 3:      # one scaling class; a few representatives are enough
        seen_special += 1
        Ds = all_hyperplanes(S)
    for D in Ds:
        if rank_gf2([abc(v)[2] for v in D]) != 8: continue
        rk = tuple(ranks_for(D, range(1, 16)))
        if rk not in cands: cands[rk] = D
print(len(cands), "distinct rank vectors; totals for n=16:", sorted(set(sum(k) for k in cands))[:6])
pickle.dump(cands, open('/tmp/linz_cands.pkl', 'wb'))
table = {}
for n in range(5, 17):
    per_node = []
    for s in range(n):
        deltas = [j ^ s for j in range(n) if j != s]
        best = min(cands.items(), key=lambda kv: sum(kv[0][d - 1] for d in deltas))
        per_node.append(sum(best[0][d - 1] for d in deltas))
        table[(n, s)] = (per_node[-1], best[1])
    print(f"n={n:2d}  per-node {per_node}  max {max(per_node)}")
pickle.dump(table, open('/tmp/linz_table.pkl', 'wb'))
