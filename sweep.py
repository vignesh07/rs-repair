"""Run searches for ISA-L point sets and keep the best verified scheme per (n, r, node)."""
import sys, re, json, subprocess, os, time
from multiprocessing import Pool
sys.path.insert(0, '/home/claude/rs')

def node_classes(n):
    m = 0
    while n % (1 << (m + 1)) == 0: m += 1
    step = 1 << m                      # A = {0..n-1} is invariant under XOR with t < 2^m
    return [s for s in range(n) if s % step == 0], step

def run_job(job):
    kind, n, r, node, iters, seed, restarts, ps = job
    pts = [str(i) for i in range(n)]
    if kind == 'sa2':
        args = ['./sa2', str(n), str(r), str(node), str(iters), str(seed), '1.0', '0.12', str(restarts), str(ps)] + pts
    else:
        args = ['./sa3', str(n), str(node), str(iters), str(seed), '1.0', '0.12', str(restarts), str(ps)] + pts
    out = subprocess.run(args, capture_output=True, text=True).stdout
    m = re.search(r'bw=(\d+)\s+(c|ab)=([\d,]+)', out)
    return {"kind": kind, "n": n, "r": r, "node": node, "bw": int(m.group(1)), "fmt": m.group(2), "vals": list(map(int, m.group(3).split(',')))}

if __name__ == '__main__':
    which = sys.argv[1]; outpath = sys.argv[2]; procs = int(sys.argv[3])
    jobs = []
    if which in ('r2', 'r2a', 'r2b'):
        lo, hi = {'r2': (4, 17), 'r2a': (4, 13), 'r2b': (13, 17)}[which]
        for n in range(lo, hi):
            for s in node_classes(n)[0]:
                jobs.append(('sa2', n, 2, s, 1500000, 100 + s, 10, 0.5))
    elif which == 'r3':
        for n in range(5, 17):
            for s in node_classes(n)[0]:
                jobs.append(('sa2', n, 3, s, 4000000, 200 + s, 8, 0.8))
    elif which == 'bb':                      # Backblaze 17+3 and a few longer three-parity codes (points 0..n-1)
        for n in (20,):
            for s in node_classes(n)[0]:
                jobs.append(('sa2', n, 3, s, 6000000, 500 + s, 4, 0.8))
    elif which == 'r4small':
        for n in range(6, 12):
            for s in node_classes(n)[0]:
                jobs.append(('sa2', n, 4, s, 5000000, 300 + s, 6, 0.8))
    elif which == 'r4lin':
        for n in range(9, 17):
            for s in node_classes(n)[0]:
                jobs.append(('sa3', n, 4, s, 8000000, 400 + s, 5, 0.8))
    t0 = time.time()
    with Pool(procs) as pool, open(outpath, 'a') as f:
        for res in pool.imap_unordered(run_job, jobs):
            f.write(json.dumps(res) + "\n"); f.flush()
            print(f"[{time.time()-t0:6.0f}s] n={res['n']} r={res['r']} node={res['node']} bw={res['bw']}", flush=True)
