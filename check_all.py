#!/usr/bin/env python3
"""Re-verify every scheme for the points 0..n-1 (schemes_isal.json, schemes_rs6_2_optimal.json, schemes_rs20_17.json).
For each scheme: check the 8 polynomials are dual codewords, recompute the bandwidth, and repair random codewords
bit-for-bit on (a) the Intel ISA-L Cauchy code and (b) the Backblaze / klauspost Vandermonde code.
Where two files hold a scheme for the same node, the cheaper one counts.
usage: python3 check_all.py [trials_per_scheme]"""
import json, sys, os
from verify import verify
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 50
best = {}; count = 0
for path in ("schemes_isal.json", "schemes_rs6_2_optimal.json", "schemes_rs20_17.json"):
    if not os.path.exists(path):
        continue
    for name, sch in json.load(open(path)).items():
        for code in ("isal", "vandermonde"):
            s = dict(sch); s["code"] = code
            bw, _ = verify(s, trials=trials, seed=count, quiet=True)
            assert bw == sch["bandwidth_bits"], (name, code, bw)
        key = (name.split(" node")[0], sch["node"])
        best[key] = min(best.get(key, 999), sch["bandwidth_bits"]); count += 1
print(f"{count} schemes verified on both codes ({trials} random repairs each)")
worst = {}
for (code, node), bw in best.items():
    worst[code] = max(worst.get(code, 0), bw)
for key, bw in worst.items():
    n, k = map(int, key[3:-1].split(","))
    print(f"{key:10s} worst node {bw:3d} bits   classical {8*k}")
