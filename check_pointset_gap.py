#!/usr/bin/env python3
"""Certify the RS(6,2) comparison in RS6_2_POINTSET_GAP.md. Standard library only."""
from collections import Counter
from functools import lru_cache
from itertools import combinations

from verify import inv, mul, power, rank_gf2, trace


@lru_cache(maxsize=None)
def divide(x, y):
    return mul(x, inv(y))


def slopes(points, target, reference):
    return {i: divide(i ^ reference, i ^ target) for i in points if i not in (target, reference)}


def reference_for(points, target):
    return next(i for i in points if i != target)


def build_protocol(points, target):
    h = reference_for(points, target)
    c = slopes(points, target, h)
    if rank_gf2(c.values()) == 4:
        beta = next(b for b in range(2, 256)
                    if rank_gf2([*c.values(), *(mul(b, x) for x in c.values())]) == 8)
        reference_basis = [1, beta]
    else:
        # This specific basis suffices for both exceptional nodes. Check rank below.
        reference_basis = [1, 2, 4]

    # A wire (i, v) carries Tr(v*f(i)). Wires from h come first. For each
    # other wire, cancellation records which reference bit removes its f(h) term.
    wires = [(h, v) for v in reference_basis]
    target_basis, cancellation = [], []
    for i, slope in c.items():
        b_i = divide(i ^ target, h ^ target)
        for ref_bit, v in enumerate(reference_basis):
            target_coeff = mul(slope, v)
            if rank_gf2(target_basis + [target_coeff]) > len(target_basis):
                target_basis.append(target_coeff)
                wires.append((i, divide(v, b_i)))
                cancellation.append(ref_bit)
    assert len(target_basis) == 8
    return wires, target_basis, cancellation


def verify_all_messages(points, target, wires, target_basis, cancellation):
    # The decoder sees only the packed wire bits. Invert the eight target traces.
    target_traces = {sum(trace(mul(v, x)) << j for j, v in enumerate(target_basis)): x
                     for x in range(256)}
    assert len(target_traces) == 256
    nref = len(wires) - 8
    decoder = []
    for transcript in range(1 << len(wires)):
        trace_bits = transcript >> nref
        for j, ref_bit in enumerate(cancellation):
            trace_bits ^= ((transcript >> ref_bit) & 1) << j
        decoder.append(target_traces[trace_bits])

    # Encode the contributions of the constant and linear terms separately.
    # For f(x)=a+b*x their XOR is the exact transcript sent by the helpers.
    constant, linear = [], []
    for x in range(256):
        constant.append(sum(trace(mul(v, x)) << j for j, (_, v) in enumerate(wires)))
        linear.append(sum(trace(mul(v, mul(x, i))) << j for j, (i, v) in enumerate(wires)))
    for b, encoded_linear in enumerate(linear):
        bt = mul(b, target)
        for a, encoded_constant in enumerate(constant):
            assert decoder[encoded_constant ^ encoded_linear] == a ^ bt


def check_profile(label, points, expected_ranks, expected_costs):
    print(f'{label}, points = {list(points)}')
    for t, rank, cost in zip(points, expected_ranks, expected_costs):
        h = reference_for(points, t)
        coefficients = list(slopes(points, t, h).values())
        assert rank_gf2(coefficients) == rank
        wires, basis, cancellation = build_protocol(points, t)
        assert len(wires) == cost
        verify_all_messages(points, t, wires, basis, cancellation)
        print(f'  t={t:3}, h={h}, slopes={coefficients}, rank={rank}, bits={cost}')


def main():
    if not __debug__:
        raise SystemExit('Run without -O so certificate assertions are enabled.')
    z = 78
    assert power(z, 4) ^ z ^ 1 == 0 and len({power(z, i) for i in range(15)}) == 15
    prescribed = (0, 1, z, power(z, 2), power(z, 3), power(z, 4))
    check_profile('ISA-L points', tuple(range(6)), [4]*6, [10]*6)
    check_profile('Prescribed subfield points', prescribed,
                  [4, 4, 4, 3, 4, 3], [10, 10, 10, 11, 10, 11])
    print('786,432 exact repairs verified (65,536 messages x 6 failures x 2 evaluation codes).')

    # For subfield points, rank 4 is necessary and sufficient for 10 bits.
    subfield = [x for x in range(256) if power(x, 16) == x]
    histogram = Counter()
    for points in combinations(subfield, 6):
        bad = sum(rank_gf2(slopes(points, t, reference_for(points, t)).values()) < 4
                  for t in points)
        histogram[bad] += 1
    assert histogram == {0: 880, 2: 5280, 4: 1320, 5: 528}
    print('Of the 8,008 subfield six-sets, 880 admit 10 bits at every node.')


if __name__ == '__main__':
    main()
