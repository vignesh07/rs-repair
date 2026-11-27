#!/usr/bin/env python3
"""Independent check of the wire format and of the traffic figures.

The C tool `rsr helper` writes the payload files.  This script then rebuilds the lost shard *from those files alone*,
in Python, using only the table in rsr_tables.h and the wire format as documented in rsrepair.h, and compares the
result with the shard the other library wrote.  It shares no code with rsrepair.c.

usage: python3 check_wire.py <stripe file> <plain|cauchy> n k        (stripe file = n shards back to back; needs ./rsr)
"""
import os, re, subprocess, sys, tempfile

path, flavour, n, k = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
HERE = os.path.dirname(os.path.abspath(__file__))


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        if a & 0x100:
            a ^= 0x11D
        b >>= 1
    return r


def inv(a):
    r = 1
    for _ in range(254):
        r = mul(r, a)
    return r


src = open(os.path.join(HERE, "rsr_tables.h")).read()
codes = [tuple(map(int, t)) for t in re.findall(r"\{(\d+), (\d+), (\d+)\}", re.search(r"rsr_codes\[.*?\{(.*?)\};", src, re.S).group(1))]
ints = lambda name: list(map(int, re.search(name + r"\[\d+\] = \{(.*?)\};", src, re.S).group(1).replace("\n", " ").strip(" ,").split(",")))
node_off, blob = ints("rsr_node_off"), ints("rsr_blob")
first = next(f for a, b, f in codes if (a, b) == (n, k))

data = open(path, "rb").read()
L = len(data) // n
shard = [data[i * L:(i + 1) * L] for i in range(n)]
def u_lost(lost):
    """chunk i of the Cauchy code holds u_i * f(i), u_i = 1 / prod_{t<k, t!=i} (i ^ t); the plain code has u = 1."""
    if flavour == "plain":
        return 1
    prod = 1
    for t in range(k):
        if t != lost:
            prod = mul(prod, lost ^ t)
    return inv(prod)

total_payload = 0
for lost in range(n):
    p = node_off[first + lost]
    fin = [mul(u_lost(lost), f) for f in blob[p:p + 8]]
    p += 8
    r, dmask = {}, {}
    for i in range(n):
        if blob[p] == 0xFF:
            p += 1
            continue
        r[i] = blob[p]
        dmask[i] = blob[p + 1 + r[i]:p + 1 + 2 * r[i]]
        p += 1 + 2 * r[i]
    with tempfile.TemporaryDirectory() as d:
        pay = {}
        for i in r:
            open(f"{d}/s", "wb").write(shard[i])
            subprocess.run([os.path.join(HERE, "rsr"), "helper", flavour, str(n), str(k), str(lost), str(i), f"{d}/s", f"{d}/p"], check=True)
            pay[i] = open(f"{d}/p", "rb").read()
            assert len(pay[i]) == r[i] * (64 * (L // 512) + (L % 512 + 7) // 8), "payload size differs from the documented formula"
            total_payload += len(pay[i])
    out = bytearray(L)
    pos = {i: 0 for i in r}                       # read offset in each payload
    for blk in range(0, L, 512):
        m = min(512, L - blk)
        pb = (m + 7) // 8
        z = [0] * m                               # intermediate byte per chunk byte
        for i in r:
            for j in range(r[i]):
                plane = pay[i][pos[i] + j * pb: pos[i] + (j + 1) * pb]
                for t in range(m):
                    if plane[t >> 3] >> (t & 7) & 1:
                        z[t] ^= dmask[i][j]
            pos[i] += r[i] * pb
        for t in range(m):
            v = 0
            for b in range(8):
                if z[t] >> b & 1:
                    v ^= fin[b]
            out[blk + t] = v
    assert bytes(out) == shard[lost], f"python rebuild differs for lost shard {lost}"
print(f"OK RS({n},{k}) {flavour}: all {n} shards rebuilt in Python from the C helper payloads; "
      f"average download {total_payload / n / (k * L):.4f} of normal repair")
