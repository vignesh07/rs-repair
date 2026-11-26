"""The searches of sweep.py with the same settings, on the F16-based points of the 2022 table / 2024 thesis:
A = {0, 1, w, w^2, ..., w^(n-2)} with w = z^17 (z = 2 in GF(256)/0x11D), exactly as printed in the thesis appendix.

Affine maps x -> beta*x + gamma that preserve A carry a scheme for one node to a scheme for another with the same
bandwidth (thesis, Proposition 6), so one node per orbit is searched and the result is transported.

usage: python3 sweep_f16.py <jobset> <out.jsonl> <procs>      jobsets: r2 r3 r3b r4small r4lin r4big all
       python3 sweep_f16.py time                              (prints the job list sizes)"""
import sys, re, json, subprocess, os, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify import mul, inv, power


def f16_points(n):
    return [0] + [power(2, 17 * e) for e in range(n - 1)]


def affine_group(pts):
    A = set(pts); maps = []
    a0, a1 = pts[0], pts[1]
    for b0 in pts:
        for b1 in pts:
            if b0 == b1:
                continue
            beta = mul(b0 ^ b1, inv(a0 ^ a1)); gamma = b0 ^ mul(beta, a0)
            if all((mul(beta, a) ^ gamma) in A for a in pts):
                maps.append((beta, gamma))
    return maps


def node_orbits(pts):
    """[(representative node, {node: beta of a map sending the representative to it})]"""
    maps = affine_group(pts); idx = {a: i for i, a in enumerate(pts)}
    seen = set(); orbits = []
    for s in range(len(pts)):
        if s in seen:
            continue
        orb = {s: 1}
        for beta, gamma in maps:
            t = idx[mul(beta, pts[s]) ^ gamma]
            orb.setdefault(t, beta)
        seen |= set(orb); orbits.append((s, orb))
    return orbits


def run_job(job):
    kind, n, r, node, iters, seed, restarts, ps = job
    pts = [str(p) for p in f16_points(n)]
    if kind == 'sa2':
        args = [os.path.join(HERE, 'sa2'), str(n), str(r), str(node), str(iters), str(seed), '1.0', '0.12', str(restarts), str(ps)] + pts
    else:
        args = [os.path.join(HERE, 'sa3'), str(n), str(node), str(iters), str(seed), '1.0', '0.12', str(restarts), str(ps)] + pts
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True).stdout
    m = re.search(r'bw=(\d+)\s+(c|ab)=([\d,]+)', out)
    return {"points": "f16", "kind": kind, "n": n, "r": r, "node": node, "seed": seed, "bw": int(m.group(1)), "fmt": m.group(2),
            "vals": list(map(int, m.group(3).split(','))), "sec": round(time.time() - t0, 1)}


def jobs_for(which):
    jobs = []
    def reps(n):
        return [s for s, _ in node_orbits(f16_points(n))]
    if which in ('r2', 'all'):
        for n in range(4, 17):
            jobs += [('sa2', n, 2, s, 1500000, 100 + s, 10, 0.5) for s in reps(n)]
    if which in ('r3', 'all'):
        for n in range(5, 17):
            jobs += [('sa2', n, 3, s, 4000000, 200 + s, 8, 0.8) for s in reps(n)]
    if which in ('r4small', 'all'):
        for n in range(6, 12):
            jobs += [('sa2', n, 4, s, 5000000, 300 + s, 6, 0.8) for s in reps(n)]
    if which in ('r4lin', 'all'):
        for n in range(9, 17):
            jobs += [('sa3', n, 4, s, 8000000, 400 + s, 5, 0.8) for s in reps(n)]
    if which in ('r3b', 'all'):
        for n in range(13, 17):
            jobs += [('sa2', n, 3, s, 4000000, 600 + s, 6, 0.8) for s in reps(n)]
    if which in ('r4big', 'all'):            # beyond the ISA-L sweep: general search for four parities, n = 12..16
        for n in range(12, 17):
            jobs += [('sa2', n, 4, s, 5000000, 700 + s, 6, 0.8) for s in reps(n)]
    return jobs


if __name__ == '__main__':
    if sys.argv[1] == 'time':
        for w in ('r2', 'r3', 'r3b', 'r4small', 'r4lin', 'r4big'):
            print(w, len(jobs_for(w)))
        for n in range(4, 17):
            print(n, [(s, len(o)) for s, o in node_orbits(f16_points(n))])
        sys.exit()
    which, outpath, procs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    done = set()
    if os.path.exists(outpath):
        for line in open(outpath):
            d = json.loads(line); done.add((d["kind"], d["n"], d["r"], d["node"], d["seed"]))
    jobs = [j for w in which.split(',') for j in jobs_for(w)]
    jobs = [j for j in jobs if (j[0], j[1], j[2], j[3], j[5]) not in done]
    print(f"{len(jobs)} jobs to run", flush=True)
    t0 = time.time()
    with Pool(procs) as pool, open(outpath, 'a') as f:
        for res in pool.imap_unordered(run_job, jobs):
            f.write(json.dumps(res) + "\n"); f.flush()
            print(f"[{time.time()-t0:6.0f}s] {res['kind']} n={res['n']} r={res['r']} node={res['node']} bw={res['bw']} ({res['sec']}s)", flush=True)
