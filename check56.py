#!/usr/bin/env python3
"""RS(16,12) over GF(256) on the 16-point subspace V = span(1, 2, 8, 20): one scheme, 56 bits for every node.
Found with the adjoint-space criterion (adj.c, mode B). Checked here with the independent verifier (real repairs)."""
import json, sys
from verify import mul, verify

ABC = [(178, 180, 1), (10, 183, 2), (104, 145, 4), (56, 202, 8), (197, 95, 16), (47, 74, 32), (190, 54, 64), (207, 63, 128)]
BASIS = [1, 2, 8, 20]
V = sorted({(a * 1) ^ (b * 2) ^ (c * 8) ^ (d * 20) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)})
assert len(V) == 16

trials = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
schemes = []
for node, t in enumerate(V):
    t2 = mul(t, t); t3 = mul(t2, t)
    # g(x) = a (x+t)^3 + b (x+t) + c, lowest power first
    polys = [[mul(a, t3) ^ mul(b, t) ^ c, mul(a, t2) ^ b, mul(a, t), a] for a, b, c in ABC]
    sch = {"n": 16, "k": 12, "points": V, "node": node, "polys": polys, "code": "eval"}
    bw, ranks = verify(sch, trials=trials, quiet=True)
    assert bw == 56, bw
    print(f"node at point {t:2d}: {bw} bits, helper bits {[ranks[j] for j in range(16) if j != node]}, {trials} repairs correct")
    schemes.append(sch)
json.dump(schemes, open("schemes_rs16_12_subspace.json", "w"))
print("all 16 nodes: 56 bits (classical 96, best published 60)")
