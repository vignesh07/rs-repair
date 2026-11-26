import sys, json, itertools
sys.path.insert(0, '/home/claude/rs')
from verify import mul, power, inv, rank_gf2, verify
from linz import subspaces_dim4, S_of, abc, ranks_for

def scheme_from_abc(S, n, node):
    """S: list of (a,b,c) giving g(y) = a y^3 + b y + c with y = x + a*, a* = node (ISA-L points). Returns verify.py scheme."""
    t = node
    polys = []
    for (a, b, c) in S:
        # expand a (x+t)^3 + b (x+t) + c  in char 2: (x+t)^3 = x^3 + t x^2 + t^2 x + t^3
        t2 = mul(t, t); t3 = mul(t2, t)
        polys.append([mul(a, t3) ^ mul(b, t) ^ c, mul(a, t2) ^ b, mul(a, t), a])
    return {"n": n, "k": n - 4, "points": list(range(n)), "node": node, "polys": polys, "code": "isal"}

if __name__ == '__main__':
    best = json.load(open('/tmp/linz_best.json'))
    S = [tuple(x) for x in best["S"]]
    for node in (0, 5, 15):
        verify(scheme_from_abc(S, 16, node), trials=500)
    # the same polynomials on the 14-point code: helpers are the 13 other points of {0..13}
    for node in range(14):
        bw, _ = verify(scheme_from_abc(S, 14, node), trials=200, quiet=True)
        print("RS(14,10) node", node, "bandwidth", bw)
