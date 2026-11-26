#!/usr/bin/env python3
"""The whole trick on one screen: rebuild a lost byte of RS(4,2) from 12 bits instead of 16.

RS(4,2) as Intel ISA-L builds it: data bytes on nodes 0 and 1, parity bytes on nodes 2 and 3.
Node 0 is lost. The classical rebuild reads two whole surviving bytes (16 bits).
Here each of the three survivors sends 4 parity bits of its byte (12 bits in total).
The three tables below depend only on the code, never on the data. They are the scheme "RS(4,2) node 0"
of schemes_isal.json, written out as bit masks."""
import random


def gf_mul(a, b):                       # product of two bytes in GF(256), polynomial 0x11D
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        if a & 0x100:
            a ^= 0x11D
        b >>= 1
    return r


def parity(x):                          # XOR of the bits of x
    return bin(x).count("1") & 1


# encode one stripe: two data bytes, two parity bytes (the constants are ISA-L's Cauchy matrix for n = 4, k = 2)
PARITY_ROWS = [[142, 244], [244, 142]]
data = [random.randrange(256), random.randrange(256)]
stripe = data + [gf_mul(row[0], data[0]) ^ gf_mul(row[1], data[1]) for row in PARITY_ROWS]

# the repair scheme for node 0
SEND = {1: [230, 215, 131, 191], 2: [245, 87, 220, 158], 3: [95, 217, 117, 71]}          # survivor j sends parity(byte & mask) for each of its 4 masks
COMBINE = [[(1, 0), (2, 0), (3, 0)], [(1, 1), (2, 1), (3, 1)], [(1, 2), (2, 2), (3, 2)], [(1, 3), (2, 3), (3, 0), (3, 2)], [(1, 0), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (3, 3)], [(1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 3), (3, 0)], [(1, 0), (1, 2), (1, 3), (2, 1), (2, 3), (3, 0), (3, 2)], [(1, 0), (2, 1), (3, 0), (3, 1), (3, 3)]]    # fact i about the lost byte = XOR of these (survivor, which of its bits)
LOST_MASKS = [127, 63, 31, 15, 135, 67, 161, 208]    # fact i is parity(lost byte & LOST_MASKS[i])

# each survivor computes 4 bits from its own byte and sends them
sent = {j: [parity(stripe[j] & m) for m in SEND[j]] for j in (1, 2, 3)}

# the rebuilder turns the 12 received bits into 8 facts about the lost byte, then solves for the byte
facts = []
for row in COMBINE:
    bit = 0
    for j, q in row:
        bit ^= sent[j][q]
    facts.append(bit)
candidates = [x for x in range(256) if all(parity(x & LOST_MASKS[i]) == facts[i] for i in range(8))]

print("stripe (nodes 0..3):", stripe, "  node 0 is lost")
for j in (1, 2, 3):
    print(f"node {j} holds {stripe[j]:3d} and sends bits {sent[j]}")
print("8 facts about the lost byte:", facts)
print("bytes that fit all 8 facts:", candidates, "  lost byte was", stripe[0])
assert candidates == [stripe[0]]
print("rebuilt from", sum(len(v) for v in sent.values()), "bits; the classical rebuild reads 16")
