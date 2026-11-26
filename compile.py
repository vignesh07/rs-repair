"""Collect the best scheme per (n, r, node) for ISA-L points, verify each end to end, and print the comparison table."""
import sys, json, glob, pickle, os
sys.path.insert(0, '/home/claude/rs')
from verify import verify, scheme_from_c
from linz import abc
from linz2 import scheme_from_abc
from sweep import node_classes

# Dinh et al., ISIT 2022 (arXiv:2205.11015), Tables II-IV: n -> (classical, ISA-L heuristic, F16 algebraic, F16 heuristic); max over nodes, bits
PUB = {
    2: {4: (16, 12, 18, 12), 5: (24, 18, 24, 18), 6: (32, 24, 30, 24), 7: (40, 32, 36, 30), 8: (48, 38, 42, 38), 9: (56, 44, 48, 44),
        10: (64, 50, 54, 50), 11: (72, 58, 60, 56), 12: (80, 64, 66, 64), 13: (88, 72, 72, 70), 14: (96, 80, 78, 76), 15: (104, 84, 84, 84), 16: (112, 90, 90, 90)},
    3: {4: (8, 8, 18, 8), 5: (16, 12, 24, 12), 6: (24, 16, 30, 16), 7: (32, 22, 36, 22), 8: (40, 28, 42, 28), 9: (48, 34, 48, 32),
        10: (56, 40, 54, 40), 11: (64, 46, 60, 46), 12: (72, 52, 66, 52), 13: (80, 58, 72, 58), 14: (88, 66, 78, 64), 15: (96, 72, 84, 70), 16: (104, 78, 90, 76)},
    4: {5: (8, 8, 16, 8), 6: (16, 12, 20, 12), 7: (24, 16, 24, 16), 8: (32, 22, 28, 22), 9: (40, 28, 32, 26), 10: (48, 36, 36, 32),
        11: (56, 42, 40, 38), 12: (64, 48, 44, 44), 13: (72, 54, 48, 48), 14: (80, 62, 52, 52), 15: (88, 68, 56, 56), 16: (96, 60, 60, 60)},
}

def load_candidates():
    best = {}      # (n, r, node) -> (bw, scheme)
    def offer(n, r, node, sch):
        key = (n, r, node)
        bw, _ = verify(sch, trials=40, quiet=True)          # cheap check now, full check for the winners later
        if key not in best or bw < best[key][0]:
            best[key] = (bw, sch)
    for path in glob.glob('/home/claude/rs/res_*.jsonl'):
        for line in open(path):
            d = json.loads(line)
            n, r, rep = d["n"], d["r"], d["node"]
            _, step = node_classes(n)
            for t in range(step):
                node = rep ^ t
                if d["fmt"] == "c":
                    sch = scheme_from_c(n, r, node, d["vals"])
                else:                                             # a_i, b_i pairs: g_{e_i} = a y^3 + b y + e_i
                    flat = []
                    for i in range(8):
                        a, b = d["vals"][2 * i], d["vals"][2 * i + 1]
                        flat += [b, 0, a]
                    sch = scheme_from_c(n, 4, node, flat)
                offer(n, r, node, sch)
    if os.path.exists('/tmp/linz_table.pkl'):
        table = pickle.load(open('/tmp/linz_table.pkl', 'rb'))
        for (n, s), (bw, D) in table.items():
            offer(n, 4, s, scheme_from_abc([abc(v) for v in D], n, s))
    return best

if __name__ == '__main__':
    best = load_candidates()
    out = {}
    rows = []
    for r in (2, 3, 4):
        for n in sorted(PUB[r]):
            keys = [(n, r, s) for s in range(n)]
            if not all(k in best for k in keys):
                continue
            per = []
            for k in keys:
                bw, sch = best[k]
                got, ranks = verify(sch, trials=1000, quiet=True)      # full end-to-end check of every reported scheme
                assert got == bw
                per.append(bw)
                out[f"RS({n},{n-r}) node {k[2]}"] = {"bandwidth_bits": bw, "helper_ranks": [ranks[j] for j in sorted(ranks)], **sch}
            cl, isal, f16a, f16h = PUB[r][n]
            pub_best = min(isal, f16h, f16a)
            flag = "NEW<isal,any" if max(per) < pub_best else ("NEW<isal" if max(per) < isal else ("=" if max(per) == isal else "worse"))
            rows.append((r, n, cl, isal, min(f16a, f16h), max(per), sum(per) / len(per), flag, per))
    print(f"{'code':10s} {'classic':>7s} {'pub ISA-L':>9s} {'pub any':>7s} {'ours max':>8s} {'ours avg':>8s}  status")
    for r, n, cl, isal, f16, mx, avg, flag, per in rows:
        print(f"RS({n},{n-r})".ljust(10), f"{cl:7d} {isal:9d} {f16:7d} {mx:8d} {avg:8.1f}  {flag}   per-node {per}")
    json.dump(out, open('/home/claude/rs/schemes_isal.json', 'w'))
    json.dump([{"r": r, "n": n, "classical": cl, "published_isal": isal, "published_best_any_points": f16, "ours_max": mx, "ours_avg": round(avg, 2), "per_node": per}
               for r, n, cl, isal, f16, mx, avg, flag, per in rows], open('/home/claude/rs/table_isal.json', 'w'), indent=1)
    print(len(out), "schemes written to schemes_isal.json (each re-verified with 1000 real repairs)")
