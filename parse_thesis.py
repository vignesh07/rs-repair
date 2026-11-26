#!/usr/bin/env python3
"""Recompute the repair bandwidth of every scheme printed in Appendix A of T. X. Dinh's thesis
(Efficient Repair of Reed-Solomon Codes, RMIT 2024) from the extracted text, and compare with this work.

usage: pdftotext -layout thesis.pdf thesis.txt && python3 -I parse_thesis.py thesis.txt table_isal.json
The thesis text is treated purely as data. z = 2 in GF(256)/0x11D (the thesis lists the ISA-L points as 0,1,z,z^25,...)."""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # verify.py sits next to this script (needed under python3 -I)
from verify import mul, power, rank_gf2, poly_eval

path = sys.argv[1]
lines = open(path, encoding="utf-8", errors="replace").read().split("\n")
start = next(i for i, l in enumerate(lines) if "Records of Lowest-Bandwidth Repair Schemes" in l and i > 6000)
end = next(i for i, l in enumerate(lines) if l.strip().startswith("Appendix B: Credits") and i > start)
body = lines[start + 1:end]

# drop page-number lines; re-attach an exponent that the PDF layout put on the line above "[x ..."
clean = []
for l in body:
    l = re.sub(r"[^\x20-\x7e∗•]", " ", l)          # strip control / private-use glyphs (pieces of the big brackets)
    s = l.strip()
    if re.fullmatch(r"\d{3}", s):
        continue
    clean.append(l)
# Three layout glitches of the text extraction are repaired here; each repair is validated below, because a wrongly
# rebuilt scheme would not have full rank at its node or would miss the stated bandwidth.
#  (1) the first polynomial of a record loses its superscripts to the line above ("2 17" above "[x +z *x+z 68"
#      means x^2 + z^17*x + z^68);
#  (2) one scheme of RS(11,8) ISA-L has a raised fragment printed on the line above;
#  (3) one scheme of RS(5,2) F16 prints "," where "*" is meant.
out = []
i = 0
while i < len(clean):
    s = clean[i].strip()
    m = re.fullmatch(r"([23])(?:\s+(\d+))?", s)
    if m:
        j = i + 1
        while j < len(clean) and not clean[j].strip():
            j += 1
        if j < len(clean) and clean[j].lstrip().startswith("[x"):
            e1, e2 = m.group(1), m.group(2)
            if e2:
                clean[j], cnt = re.subn(r"\[x\s*\+?\s*z\s*∗", f"[x{e1} + z {e2} ∗", clean[j], count=1)
                assert cnt == 1, clean[j][:80]
            else:
                clean[j] = clean[j].replace("[x", "[x" + e1, 1)
            i += 1
            continue
    if "201 ∗x+z 54" in clean[i] and i + 1 < len(clean) and "z 31,z ∗x +z" in clean[i + 1]:
        clean[i + 1], cnt = re.subn(r"z 31,z ∗x \+z\s*,", "z 31 , z 3 ∗ x2 + z 201 ∗ x + z 54 , z 3 ∗ x2 + z 3 ∗ x + z 227 , z 85 ∗ x2 + z 85 ∗ x ,", clean[i + 1])
        assert cnt == 1
        i += 1
        continue
    if "z 35 ∗ x2 + z 18 , x + z 86" in clean[i]:
        clean[i] = clean[i].replace("z 35 ∗ x2 + z 18 , x + z 86", "z 35 ∗ x2 + z 18 ∗ x + z 86")
    out.append(clean[i])
    i += 1
text = " ".join(out)
text = re.sub(r"[^\x20-\x7e∗•]", " ", text)          # strip the private-use glyphs of the big brackets

records = re.split(r"•\s*RS\(", text)[1:]


def parse_elem(tok):
    tok = tok.strip()
    if tok == "0":
        return 0
    if tok == "1":
        return 1
    m = re.fullmatch(r"z\s*(\d+)?", tok)
    if not m:
        raise ValueError(f"bad element {tok!r}")
    return power(2, int(m.group(1) or 1))


TERM = re.compile(r"(?:(z)\s*(\d+)?|(1))?\s*(?:∗)?\s*(?:(x)\s*(\d)?)?")


def parse_poly(s, r):
    coeffs = [0] * r
    for term in s.split("+"):
        t = term.strip()
        m = TERM.fullmatch(t)
        if not m or not t:
            raise ValueError(f"bad term {t!r} in {s!r}")
        c = power(2, int(m.group(2) or 1)) if m.group(1) else 1
        d = int(m.group(5) or 1) if m.group(4) else 0
        if d >= r:
            raise ValueError(f"degree {d} too high in {s!r}")
        coeffs[d] ^= c
    return coeffs


results = {}
problems = []
nolist = []
for rec in records:
    m = re.match(r"\s*(\d+),\s*(\d+)\),\s*(ISA-L|F16\s*-based)", rec)
    n, k = int(m.group(1)), int(m.group(2))
    kind = "isal" if m.group(3).startswith("ISA") else "f16"
    r = n - k
    lists = re.findall(r"\[([^\[\]]*)\]", rec)
    if "x" in lists[0]:                      # some records print no point list: use the definition from Section 5.1.3
        pts = list(range(n)) if kind == "isal" else [0] + [power(2, 17 * e) for e in range(n - 1)]
        nolist.append((n, k, kind))
    else:
        pts = [parse_elem(t) for t in lists[0].split(",")]
        lists = lists[1:]
    stated = re.search(r"repair bandwidth is\s*(\d+)", rec)
    stated = int(stated.group(1)) if stated else None
    schemes = lists
    key = (n, k, kind)
    if len(pts) != n or len(set(pts)) != n:
        problems.append((key, f"points parse: {len(pts)} points")); continue
    if kind == "isal" and pts != list(range(n)):
        problems.append((key, f"ISA-L points are not 0..n-1: {pts}"))
    prof = {}
    for idx, sch in enumerate(schemes):
        try:
            polys = [parse_poly(p, r) for p in sch.split(",")]
        except ValueError as e:
            problems.append((key, f"scheme {idx}: {e}")); continue
        if len(polys) != 8:
            problems.append((key, f"scheme {idx}: {len(polys)} polynomials")); continue
        ranks = [rank_gf2([poly_eval(g, a) for g in polys]) for a in pts]
        full = [j for j in range(n) if ranks[j] == 8]
        for j in full:                       # the scheme repairs every node where the rank is 8
            bw = sum(ranks) - 8
            if j not in prof or bw < prof[j]:
                prof[j] = bw
        if not full:
            problems.append((key, f"scheme {idx}: no node has full rank (ranks {ranks})"))
    missing = [j for j in range(n) if j not in prof]
    if missing:
        problems.append((key, f"no scheme repairs nodes {missing}"))
    results[key] = {"stated": stated, "nschemes": len(schemes), "profile": [prof.get(j) for j in range(n)],
                    "max": max(prof.values()) if prof and not missing else None, "points": pts}

ours = {}
if len(sys.argv) > 2:
    for row in json.load(open(sys.argv[2])):
        ours[(row["n"], row["n"] - row["r"])] = row
print(f"{len(records)} records parsed, {len(problems)} problems, {len(nolist)} records without a printed point list: {nolist}")
for p in problems:
    print("  PROBLEM", p)
print()
print("code       | thesis ISA-L: stated / recomputed | thesis F16: stated / recomputed | ours (ISA-L points) | verdict")
codes = sorted({(n, k) for n, k, _ in results}, key=lambda t: (t[0] - t[1], t[0]))
summary = {"record": 0, "isal_only": 0, "tie": 0, "worse": 0}
for n, k in codes:
    a = results.get((n, k, "isal")); b = results.get((n, k, "f16"))
    o = ours.get((n, k))
    ov = (o.get("ours_max") if o else None)
    best_recomputed = min(x for x in [a and a["max"], b and b["max"]] if x)
    verdict = ""
    if ov is not None:
        if ov < best_recomputed: verdict = "ours lower than every thesis scheme"; summary["record"] += 1
        elif a and a["max"] and ov < a["max"]: verdict = "ours lower on ISA-L points only"; summary["isal_only"] += 1
        elif ov == best_recomputed: verdict = "tie"; summary["tie"] += 1
        else: verdict = "THESIS LOWER"; summary["worse"] += 1
    print(f"RS({n},{k})".ljust(10), "|", f"{a['stated']} / {a['max']}".ljust(33) if a else "-".ljust(33), "|",
          f"{b['stated']} / {b['max']}".ljust(31) if b else "-".ljust(31), "|", str(ov).ljust(19), "|", verdict)
print()
print(summary)
bad = 0
for (n, k), o in ours.items():
    a = results.get((n, k, "isal")); b = results.get((n, k, "f16"))
    if not a or not b:
        print("not in thesis:", (n, k)); bad += 1; continue
    if o["published_isal"] != a["max"] or o["published_best_any_points"] != min(a["max"], b["max"]):
        print("published columns differ from thesis for", (n, k), o["published_isal"], o["published_best_any_points"], a["max"], b["max"]); bad += 1
print("published columns of our table equal the recomputed thesis values" if not bad else f"{bad} rows differ")
json.dump({f"{n},{k},{kind}": v for (n, k, kind), v in results.items()}, open(os.path.join(HERE, "thesis_appendix_recomputed.json"), "w"))
