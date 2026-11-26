"""Collect the best scheme per (n, r, node) on the F16-based points, verify each with real repairs (evaluation code on
those points), and compare with the ISA-L points and with the schemes printed in the 2024 thesis.

usage: python3 compile_f16.py res_f16.jsonl [trials]      needs table_isal.json and thesis_appendix_recomputed.json alongside"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify import verify, scheme_from_c, mul, inv, power
from sweep_f16 import f16_points, node_orbits


def to_c(d):
    if d["fmt"] == "c":
        return list(d["vals"])
    flat = []                                   # a_i, b_i pairs: g_{e_i} = a y^3 + b y + e_i
    for i in range(8):
        a, b = d["vals"][2 * i], d["vals"][2 * i + 1]
        flat += [b, 0, a]
    return flat


def transport(cflat, r, beta):
    """scheme for node a* -> scheme for node beta*a* + gamma: the coefficient of (x - a*)^m is multiplied by beta^(-m)"""
    bi = inv(beta); out = []; it = iter(cflat)
    for _ in range(8):
        for m in range(1, r):
            out.append(mul(next(it), power(bi, m)))
    return out


def main():
    path = sys.argv[1]; trials = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    orbits = {n: dict(node_orbits(f16_points(n))) for n in range(4, 17)}
    best = {}
    nlines = 0
    for line in open(path):
        d = json.loads(line); nlines += 1
        n, r, rep = d["n"], d["r"], d["node"]
        pts = f16_points(n); c = to_c(d)
        for node, beta in orbits[n][rep].items():
            sch = scheme_from_c(n, r, node, transport(c, r, beta), points=pts); sch["code"] = "eval"
            key = (n, r, node)
            if key in best and best[key][0] <= d["bw"]:
                continue
            bw, _ = verify(sch, trials=5, quiet=True)
            assert bw == d["bw"], (key, bw, d["bw"])          # the search's own count, and transport, both confirmed
            best[key] = (bw, sch)
    thesis = json.load(open(os.path.join(HERE, "thesis_appendix_recomputed.json")))
    isal = {(row["n"], row["r"]): row for row in json.load(open(os.path.join(HERE, "table_isal.json")))}
    rows = []; schemes = []
    for r in (2, 3, 4):
        for n in range(r + 2, 17):
            k = n - r
            keys = [(n, r, s) for s in range(n)]
            if not all(key in best for key in keys):
                print(f"RS({n},{k}): incomplete ({sum(key in best for key in keys)}/{n} nodes)"); continue
            per = []
            for key in keys:
                bw, sch = best[key]
                got, ranks = verify(sch, trials=trials, quiet=True)
                assert got == bw
                per.append(bw)
                schemes.append({**sch, "bandwidth": bw, "helper_ranks": [ranks[j] for j in range(n) if j != key[2]]})
            tf = thesis[f"{n},{k},f16"]; ti = thesis[f"{n},{k},isal"]; oi = isal[(n, r)]
            known_f16 = [min(a, b) for a, b in zip(per, tf["profile"])]
            known_isal = [min(a, b) for a, b in zip(oi["per_node"], ti["profile"])]
            rows.append({"r": r, "n": n, "k": k, "classical": 8 * k,
                         "thesis_isal": ti["max"], "thesis_f16": tf["max"],
                         "ours_isal": oi["ours_max"], "ours_f16": max(per), "ours_f16_avg": round(sum(per) / n, 2), "ours_f16_per_node": per,
                         "best_known_isal": max(known_isal), "best_known_f16": max(known_f16),
                         "best_known_isal_avg": round(sum(known_isal) / n, 2), "best_known_f16_avg": round(sum(known_f16) / n, 2)})
    json.dump(schemes, open(os.path.join(HERE, "schemes_f16.json"), "w"))
    json.dump(rows, open(os.path.join(HERE, "table_f16.json"), "w"), indent=1)
    print(f"{nlines} search results, {len(schemes)} schemes verified with {trials} real repairs each")
    print("code       classical | thesis: ISA-L F16 | this work: ISA-L  F16 | best known: ISA-L  F16 | lower")
    tally = {"ISA-L": 0, "F16": 0, "equal": 0}
    for w in rows:
        a, b = w["best_known_isal"], w["best_known_f16"]
        low = "ISA-L" if a < b else "F16" if b < a else "equal"
        tally[low] += 1
        print(f"RS({w['n']},{w['k']})".ljust(10), str(w["classical"]).rjust(8), "|", str(w["thesis_isal"]).rjust(12), str(w["thesis_f16"]).rjust(4), "|",
              str(w["ours_isal"]).rjust(16), str(w["ours_f16"]).rjust(4), "|", str(a).rjust(16), str(b).rjust(4), "|", low)
    print(tally)


if __name__ == "__main__":
    main()
