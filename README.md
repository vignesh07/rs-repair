# Lower repair bandwidth for deployed Reed–Solomon codes

Single-node repair schemes (Guruswami–Wootters trace repair, base field GF(2)) for Reed–Solomon codes over GF(256)
on the evaluation points **0, 1, …, n−1** (polynomial 0x11D). These are the points of Intel ISA-L's Cauchy code and of
the Backblaze / klauspost Vandermonde code, so every scheme here applies to data already encoded with those libraries:
only the repair procedure changes.

**Two things are here.** [`lib/`](lib/) is a working C library that rebuilds a lost chunk with 25–44% less download
than normal repair, at the speed of Intel ISA-L, tested on stripes written by ISA-L, klauspost/reedsolomon, Backblaze
JavaReedSolomon and Hadoop. The rest of this page is the repair schemes and proofs behind it.

Status: computational results, 26 Nov 2026. Unrefereed. Every scheme is machine-checked by `check_all.py`.

## Summary

1. **A counterexample, with proof.** Dinh et al. (ISIT 2022, Section IV-B) ask whether their prescribed GF(16)-based
   evaluation points always repair at least as cheaply as the ISA-L points. For RS(6,2) they do not: the optimal
   worst-node cost is 10 bits on the ISA-L points and 11 bits on theirs. Proof and constructions:
   [RS6_2_POINTSET_GAP.md](RS6_2_POINTSET_GAP.md); checks: `check_pointset_gap.py` and, independently, `review_gap.py`.
2. **Lower best-known repair costs** for 27 of 36 short codes on the ISA-L points, e.g. RS(14,10) at 50 bits per lost
   byte (published best 52, classical 80). Table below.
3. **Three codes proven optimal**: RS(4,2), RS(5,2), RS(6,2).
4. **The two point sets under one search.** The same programs and settings on both: the ISA-L points come out lower
   for 26 of 36 codes, the GF(16)-based points for 1, and 9 are equal. Only the RS(6,2) row is a proof; the rest
   compares best-known schemes. Section "Both point sets under the same search".
5. **Why some helpers need 3 bits instead of 4** in the four-parity subspace construction, with an exact criterion and
   a 56-bit RS(16,12) that it produces.
6. **A working library**: [`lib/`](lib/) turns the schemes for the points 0..n−1 into code. Rebuilding a
   chunk downloads 25–44% less than normal repair. The rebuild is 1.1–2.2× faster than Intel ISA-L's own
   reconstruction on chunks of 1 MiB and up and within 0.7–1.8× of it on smaller buffers, and is tested byte for byte
   on stripes written by ISA-L, klauspost/reedsolomon, Backblaze JavaReedSolomon and Hadoop's RS coder. Costs and limits are in [`lib/README.md`](lib/README.md).

## Results (bits downloaded to rebuild one lost byte; worst node unless stated)

"Published" = Dinh, Nguyen, Mohan, Boztas, Luong, Dau, *Practical Considerations in Repairing Reed-Solomon Codes*,
ISIT 2022 (arXiv:2205.11015), Tables II–IV: best value for these points, and best value for any evaluation points.

| code | classical | published, these points | published, any points | **this work** | avg over nodes | status vs. 2022 table |
|---|---|---|---|---|---|---|
| RS(4,2) | 16 | 12 | 12 | **12** | 12.0 | ties; optimal |
| RS(5,3) | 24 | 18 | 18 | **18** | 17.2 | ties |
| RS(6,4) | 32 | 24 | 24 | **24** | 24.0 | ties |
| RS(7,5) | 40 | 32 | 30 | **30** | 30.0 | new for these points |
| RS(8,6) | 48 | 38 | 38 | **35** | 35.0 | **new record** |
| RS(9,7) | 56 | 44 | 44 | **42** | 41.4 | **new record** |
| RS(10,8) | 64 | 50 | 50 | **49** | 48.0 | **new record** |
| RS(11,9) | 72 | 58 | 56 | **55** | 53.7 | **new record** |
| RS(12,10) | 80 | 64 | 64 | **60** | 59.3 | **new record** |
| RS(13,11) | 88 | 72 | 70 | **67** | 65.8 | **new record** |
| RS(14,12) | 96 | 80 | 76 | **73** | 72.3 | **new record** |
| RS(15,13) | 104 | 84 | 84 | **80** | 78.4 | **new record** |
| RS(16,14) | 112 | 90 | 90 | **84** | 84.0 | **new record** |
| RS(5,2) | 16 | 12 | 12 | **11** | 11.0 | **new record**; optimal |
| RS(6,3) | 24 | 16 | 16 | **16** | 16.0 | ties |
| RS(7,4) | 32 | 22 | 22 | **22** | 21.3 | ties |
| RS(8,5) | 40 | 28 | 28 | **27** | 27.0 | **new record** |
| RS(9,6) | 48 | 34 | 32 | **33** | 32.1 | new for these points |
| RS(10,7) | 56 | 40 | 40 | **39** | 37.6 | **new record** |
| RS(11,8) | 64 | 46 | 46 | **45** | 43.7 | **new record** |
| RS(12,9) | 72 | 52 | 52 | **50** | 49.3 | **new record** |
| RS(13,10) | 80 | 58 | 58 | **56** | 55.6 | **new record** |
| RS(14,11) | 88 | 66 | 64 | **63** | 62.0 | **new record** |
| RS(15,12) | 96 | 72 | 70 | **69** | 68.4 | **new record** |
| RS(16,13) | 104 | 78 | 76 | **75** | 75.0 | **new record** |
| RS(6,2) | 16 | 12 | 12 | **10** | 10.0 | **new record**; optimal |
| RS(7,3) | 24 | 16 | 16 | **16** | 15.1 | ties |
| RS(8,4) | 32 | 22 | 22 | **18** | 18.0 | **new record** |
| RS(9,5) | 40 | 28 | 26 | **25** | 24.2 | **new record** |
| RS(10,6) | 48 | 36 | 32 | **31** | 30.8 | **new record** |
| RS(11,7) | 56 | 42 | 38 | **38** | 36.8 | new for these points |
| RS(12,8) | 64 | 48 | 44 | **42** | 41.7 | **new record** |
| RS(13,9) | 72 | 54 | 48 | **46** | 45.5 | **new record** |
| RS(14,10) | 80 | 62 | 52 | **50** | 49.4 | **new record** |
| RS(15,11) | 88 | 68 | 56 | **54** | 53.2 | **new record** |
| RS(16,12) | 96 | 60 | 60 | **57** | 57.0 | **new record** |

Not in the published table: **RS(20,17)** (Backblaze Vaults parameters, points 0..19): 104 bits worst node
(93–104 per node) against 136 classical. This one had only a light search.

36 table rows computed: 27 beat every published value, 3 beat the published value for these points, 6 tie.

**Three rows are proven optimal.** Any repair of one node of an (n, k) MDS code over GF(256) downloads at least
8(n−1)/(n−k) bits per lost byte (the cut-set bound of regenerating codes, applied to n−1 helpers; fewer helpers only
raise it). RS(4,2) at 12 bits and RS(6,2) at 10 bits meet it exactly. RS(5,2) at 11 bits meets it after rounding up from
10.67, so it is optimal among schemes that send whole bits per byte. For RS(6,2) the 10-bit schemes for all six nodes
are in `schemes_rs6_2_optimal.json`; the search programs here had reached 10 on four nodes and 11 on two.

## Verify

    python3 check_all.py 200      # ~8 min; needs only the Python standard library
    python3 -B check_pointset_gap.py              # the RS(6,2) counterexample, under a second
    python3 review_gap.py                         # independent check of it, about a minute
    python3 compile_f16.py res_f16.jsonl 200      # the 377 schemes on the GF(16)-based points and the comparison
    python3 check56.py                            # the 56-bit RS(16,12)
    python3 toy_repair.py                         # one byte of RS(4,2) rebuilt from 12 bits

For each scheme (403 in the three files for these points) it (1) builds the generator matrix exactly as ISA-L (`gf_gen_cauchy1_matrix`) and as
Backblaze/klauspost (systematic Vandermonde) do, (2) checks the eight check polynomials are dual codewords,
(3) recomputes the bandwidth as the sum of GF(2)-ranks at the helpers, and (4) repairs random codewords using only those bits.
`verify.py` shares no code with the search programs.

Scheme format (`schemes_isal.json`, `schemes_rs20_17.json`): eight polynomials of degree < r, coefficients lowest power
first, bytes in GF(256)/0x11D. Helper j sends `helper_ranks[j]` parity bits of its byte.

## How the schemes were found

1. **General search** (`sa2.c`). A scheme for node a* is an 8-dimensional GF(2)-space D of check polynomials on which
   evaluation at a* is injective; helper j pays 8 − dim(D ∩ {g : g(a_j) = 0}). The search anneals over *all* such D, with
   proposals that force one direction of D to be a completely split polynomial. The 2022 search assumed every basis
   polynomial splits completely; dropping that assumption is where most small-code gains come from.
2. **Subspace family for four parities** (`linz.py`, `linz3.py`). The points 0..15 form a GF(2)-subspace V. For
   g(y) = a·y³ + b·y + c (y = x − a*), g(d)·d is GF(2)-linear in d. Requiring a·d⁴ + b·d² + c·d ∈ Y for all d ∈ V, with Y a
   4-dimensional subspace, gives ≤ 4 bits per helper whenever evaluation at a* stays injective (true for about a third of
   all Y on these points). This is in the style of Berman, Buzaglo, Dor, Shany, Tamo (arXiv:2012.11166). Their Theorem 17
   guarantees 4 bits per helper only when the subspace is GF(16)-linear (5 bits for a general subspace such as this one),
   but the 2022 table already lists 60 bits for these points at n = 16. Shortening such a scheme is also prior art
   (their Section VI shortens [16,12] on GF(16) to [14,10] at 52 bits; US patent 10686471, Remark 2.5), so 4(n−1) bits
   for n < 16 on these points is a corollary, not a record. What appears to be new: for particular Y three helpers need only **3** bits, giving 57 at n = 16 and
   42 / 46 / 50 / 54 at n = 12..15, each below 4(n−1).
3. **Restricted search** (`sa3.c`) inside the family g = a·y³ + b·y + c, used for n = 9..11.

## Why some helpers need only 3 bits (`adj.c`, `check56.py`)

Notation: F = GF(256), V a 4-dimensional GF(2)-subspace of points, Y a 4-dimensional subspace with trace-dual Y⊥,
S_Y = {(a,b,c) : a·d⁴ + b·d² + c·d ∈ Y for all d ∈ V}. For a subspace W write A(W) for the span of the products
∏(H∖0) over the hyperplanes H of W (the maximal Moore minors of a basis of W; dim A(W) = dim W).
For a helper at offset d put W_d = {v² + d·v : v ∈ V}, a 3-dimensional space.

**Criterion.** dim S_Y = 8 + δ₀ and the helper at offset d sends 4 + δ₀ − δ_d bits, where
δ₀ = dim{λ ∈ F : λ·A(V) ⊆ Y⊥} and δ_d = dim{λ ∈ F : λ·A(W_d) ⊆ Y⊥}.
So in the usual case dim S_Y = 8, a helper sends fewer than 4 bits exactly when Y⊥ contains a scalar multiple of A(W_d)
(3 bits in every case found; δ_d = 2 never occurred).

Sketch of proof. The members of S_Y that vanish at d are M(x² + d·x) with M(z) = αz² + βz, and the condition on them
is M(W_d) ⊆ Y: twelve GF(2)-linear conditions on (α, β). A dependency among the twelve is a triple
(y₁,y₂,y₃) ∈ (Y⊥)³ with Σ yᵢwᵢ = Σ yᵢwᵢ² = 0 for a basis w of W_d, and over F those two equations have a
one-dimensional solution space spanned by the 2×2 Moore minors of w. The same count with V and three equations gives dim S_Y.

Checked by machine: `./adj A` compares the criterion with direct computation for all 200,787 subspaces Y on the points
0..15 (3,212,592 predictions, 0 mismatches) and reproduces the counts found independently by `linz.py`.

**Consequence.** `./adj B` runs through every 16-point subspace V (up to scaling: the 11,811 containing 1) and every Y
that the criterion allows to go below 60 bits, including the cases dim S_Y = 9, 10. Best value in this family:
**56 bits** (four 3-bit helpers), reached by 4,260 of those subspaces; 5,270 reach 57 (the points 0..15 are among
them) and 2,280 reach 58. No 16-point subspace does better inside this family.
On the subfield V = GF(16) (the points used in Berman et al.'s example and in the F16-based codes of the 2022 table)
the effect cannot occur: of the 189,040 valid Y none has a 3-bit helper, so 60 bits is the best this family can do there
(`./adj A 1 152 78 10`; the criterion again matches direct computation on all 3,212,592 predictions).
Example: RS(16,12) on V = span(1, 2, 8, 20), 56 bits for every node, `python3 check56.py` (16,000 real repairs).
This is a different code from ISA-L's, so unlike the rest of this bundle it does not apply to data already encoded.

## Both point sets under the same search (`sweep_f16.py`, `compile_f16.py`)

The 2022 paper defines two families of evaluation points: the ISA-L points 0..n−1 and the GF(16)-based points
0, 1, w, …, w^(n−2), with w a root of x⁴ + x + 1 (z¹⁷ in the thesis appendix). The searches behind the ISA-L results
were rerun on the GF(16)-based points with the same programs and settings (`sa2` for two and three parities and for four
parities up to n = 11, `sa3` for four parities with n = 9..16), plus a general `sa2` run for four parities with
n = 12..16 that the ISA-L side did not get. 355 search jobs; 377 schemes, one per node, each verified with 500 real
repairs on the evaluation code at those points (`schemes_f16.json`).

"Best known" is the worst node after taking, for every node, the cheaper of the thesis's printed scheme and the scheme
found here. Bits per lost byte.

| code | classical | thesis, ISA-L points | thesis, F16 points | this search, ISA-L points | this search, F16 points | best known, ISA-L | best known, F16 | lower |
|---|---|---|---|---|---|---|---|---|
| RS(4,2) | 16 | 12 | 12 | 12 | 12 | **12** | **12** | equal |
| RS(5,3) | 24 | 18 | 18 | 18 | 18 | **18** | **18** | equal |
| RS(6,4) | 32 | 24 | 24 | 24 | 24 | **24** | **24** | equal |
| RS(7,5) | 40 | 32 | 30 | 30 | 30 | **30** | **30** | equal |
| RS(8,6) | 48 | 38 | 38 | 35 | 37 | **35** | **37** | ISA-L |
| RS(9,7) | 56 | 44 | 44 | 42 | 44 | **42** | **44** | ISA-L |
| RS(10,8) | 64 | 50 | 50 | 49 | 50 | **49** | **50** | ISA-L |
| RS(11,9) | 72 | 58 | 56 | 55 | 56 | **55** | **56** | ISA-L |
| RS(12,10) | 80 | 64 | 64 | 60 | 64 | **60** | **64** | ISA-L |
| RS(13,11) | 88 | 72 | 70 | 67 | 70 | **67** | **70** | ISA-L |
| RS(14,12) | 96 | 80 | 76 | 73 | 76 | **73** | **76** | ISA-L |
| RS(15,13) | 104 | 84 | 84 | 80 | 84 | **80** | **84** | ISA-L |
| RS(16,14) | 112 | 90 | 90 | 84 | 91 | **84** | **90** | ISA-L |
| RS(5,2) | 16 | 12 | 12 | 11 | 11 | **11** | **11** | equal |
| RS(6,3) | 24 | 16 | 16 | 16 | 16 | **16** | **16** | equal |
| RS(7,4) | 32 | 22 | 22 | 22 | 22 | **22** | **22** | equal |
| RS(8,5) | 40 | 28 | 28 | 27 | 28 | **27** | **28** | ISA-L |
| RS(9,6) | 48 | 34 | 32 | 33 | 33 | **33** | **32** | F16 |
| RS(10,7) | 56 | 40 | 40 | 39 | 40 | **39** | **40** | ISA-L |
| RS(11,8) | 64 | 46 | 46 | 45 | 46 | **45** | **46** | ISA-L |
| RS(12,9) | 72 | 52 | 52 | 50 | 52 | **50** | **52** | ISA-L |
| RS(13,10) | 80 | 58 | 58 | 56 | 57 | **56** | **57** | ISA-L |
| RS(14,11) | 88 | 66 | 64 | 63 | 64 | **63** | **64** | ISA-L |
| RS(15,12) | 96 | 72 | 70 | 69 | 70 | **69** | **70** | ISA-L |
| RS(16,13) | 104 | 78 | 76 | 75 | 76 | **75** | **76** | ISA-L |
| RS(6,2) | 16 | 12 | 12 | 10 | 11 | **10** | **11** | ISA-L |
| RS(7,3) | 24 | 16 | 16 | 16 | 16 | **16** | **16** | equal |
| RS(8,4) | 32 | 22 | 22 | 18 | 22 | **18** | **22** | ISA-L |
| RS(9,5) | 40 | 28 | 26 | 25 | 27 | **25** | **26** | ISA-L |
| RS(10,6) | 48 | 36 | 32 | 31 | 32 | **31** | **32** | ISA-L |
| RS(11,7) | 56 | 42 | 38 | 38 | 38 | **38** | **38** | equal |
| RS(12,8) | 64 | 48 | 44 | 42 | 44 | **42** | **44** | ISA-L |
| RS(13,9) | 72 | 54 | 48 | 46 | 48 | **46** | **48** | ISA-L |
| RS(14,10) | 80 | 62 | 52 | 50 | 52 | **50** | **52** | ISA-L |
| RS(15,11) | 88 | 68 | 56 | 54 | 56 | **54** | **56** | ISA-L |
| RS(16,12) | 96 | 60 | 60 | 57 | 60 | **57** | **60** | ISA-L |

ISA-L points lower: 26 codes. GF(16)-based points lower: 1 code, RS(9,6). Equal: 9 codes.

How to read this:
* RS(6,2) is proved (10 against 11, see the Summary). RS(4,2) and RS(5,2) are proved equal, both at the cut-set bound.
  Every other row compares best-known schemes and could change with a better search on either side.
* On the GF(16)-based points this search mostly lands on the thesis's own values: it is lower for 4 codes (RS(8,6) 37,
  RS(5,2) 11, RS(6,2) 11, RS(13,10) 57) and higher for 3 (RS(16,14), RS(9,6), RS(9,5)), where the thesis's schemes are kept.
* For four parities and n ≥ 12 the ISA-L values come from the 3-bit helpers of the subspace family. On the subfield
  those cannot occur (previous section), and neither search program found anything below 4(n−1) there.

## Literature check (26 Nov 2026)

A separate search of 2022–2026 work (citing papers of the 2022 table and of Berman et al., Semantic Scholar / OpenAlex
sweeps, patents, public repositories of AI-found results) found **no public scheme below 52 bits for RS(14,10) or below
60 bits for RS(16,12)** over GF(256), for any evaluation points. Closest later result found: a formula in
Liu–Xu–Zhang (arXiv:2412.18430) appears to imply 88 bits for an RS(16,14) on a special subspace (not checked
independently; 84 here, on the ISA-L points).
The public scheme files of the 2022 authors use only 2, 4 or 6 bits per helper; most schemes here use odd counts.

Nearest prior work to the 3-bit criterion: Dinh, Le, Dau, Boztas, Kruglik, Kiah, Viterbo, Etzion, Chee, *Repairing a
Single Erasure in Reed-Solomon Codes with Side Information* (ISIT 2024, arXiv:2405.07180). For the kernel-type
subspace-polynomial schemes g(x) = L_W(β(x − a*))/(x − a*) it counts each helper's bits as a dimension minus the dimension
of a subspace intersection, and optimises over the subspaces. In that construction the choice of subspace changes nothing
unless side information about the lost symbol is available, and the paper gives no numbers for short codes. The criterion
above is the same kind of count for the image-type construction used here, where the choice does matter without side
information. So the criterion extends a known style of analysis; it is not a new kind of argument.
(Checked against Section 6.1.5 of the thesis below, whose Lemma 20 states this count.)

**Thesis check.** T. X. Dinh's PhD thesis (*Efficient Repair of Reed-Solomon Codes*, RMIT, 2024,
DOI 10.25439/rmt.28192076) was read on 26 Nov 2026. Its Tables 5.1–5.3 carry the same values as the 2022 tables, and
its Appendix A prints explicit schemes for all 76 entries (38 codes, each on the ISA-L points and on F16-based points).
`parse_thesis.py` recomputes the bandwidth of every one of those schemes from the thesis text; all agree with the tables
(two appendix entries state 12 bits by a slip; their schemes give the tabulated 38 and 90). Against the recomputed
values, this work is lower than every thesis scheme for 27 codes, lower on the ISA-L points only for 3, equal for 6 and
higher for none. To repeat: `pdftotext -layout thesis.pdf thesis.txt && python3 -I parse_thesis.py thesis.txt table_isal.json`.

The thesis points to the authors' leaderboard (dausonhoang.github.io/rsbandwidth) for later updates. That page returned
errors on 26 Nov 2026, so anything added after the thesis could not be checked.

**Counterexample to the fixed-point-family comparison.** Section IV-B of the 2022 paper and Section 5.1.3 of the thesis
ask whether their prescribed F16-based points always repair at least as cheaply as the ISA-L points. For RS(6,2) the
optimal worst-node costs are **11 bits on their F16-based points and 10 bits on the ISA-L points**. A binary rank
obstruction excludes 10 bits at two of the six subfield points, and explicit protocols attain each optimum.
See the [proof and exact profiles](RS6_2_POINTSET_GAP.md) and run `python3 -B check_pointset_gap.py` to check all
786,432 message, target, and point-set combinations. Other six-point subsets of GF(16) can attain 10 bits at every node.

## Honest caveats

* "New record" means: below every published value found by the search above. It is not a proof that nothing better exists.
* Only RS(4,2), RS(5,2) and RS(6,2) are proven optimal (see above). For the rest the proven floors are far below the
  values here: for RS(14,10) the cut-set bound is 26 bits and the counting bound for linear schemes 28; for RS(16,12)
  they are 30 and 36.
* The RS(6,2) lower bound is for repair of one byte at a time, each helper sending a fixed number of bits computed from
  its own byte, which is the setting of the 2022 paper. The proof in the note uses nothing about the size of the
  symbol, so it also excludes an average of 10 bits per byte for XOR-linear repair of several bytes jointly at the two
  obstructed nodes; how close to 10 such schemes can get is open, and nonlinear joint repair is not covered. Joint
  repair does matter in general: for [4,2] MDS codes see H. Sun, *Capacity Meets the Four Subspace Problem*
  (arXiv:2610.42023, appendix on repair; read here through an automated summary).
* Bandwidth only: every helper still reads its whole byte, and all n−1 survivors take part (classical repair needs k).
* Compute cost: the thesis reports a C implementation of trace repair on top of ISA-L that is 1.8–2.8 times slower than
  ISA-L's own repair (its Table 5.4), argues that the network saving can outweigh this, and leaves a test inside a
  storage system (Hadoop) as future work. The library in [`lib/`](lib/) closes that gap for the rebuild step on the one
  machine it was measured on (its README says where ISA-L is still faster, and what the n−1 helpers spend); a test
  inside a storage system is still future work.
* The trace-repair framework (Guruswami–Wootters), the subspace idea and shortening are prior work.
* The 3-bit criterion is a short duality argument, given above as a sketch and confirmed by exhaustive computation for
  GF(256), four parities, 16 points. It has not been refereed. The literature check looked for schemes below 4 bits per
  helper and found none; it did not look for this criterion under other names (adjoint subspaces, Gabidulin duals).

## Files

`verify.py` independent verifier · `check_all.py` verifies everything on the points 0..n−1 · `toy_repair.py` the whole
trick on one screen (RS(4,2), 12 bits instead of 16) · `schemes_*.json` the schemes · `table_isal.json`
the comparison table · `sa.c`, `sa2.c`, `sa3.c`, `gf.h` search programs · `linz*.py` subspace family · `sweep.py`, `compile.py` drivers ·
`adj.c` 3-bit criterion and scan over point subspaces (`gcc -O2 -o adj adj.c && ./adj A && ./adj B`) · `check56.py`,
`schemes_rs16_12_subspace.json` the 56-bit RS(16,12) · `RS6_2_POINTSET_GAP.md`, `check_pointset_gap.py`,
`schemes_rs6_2_optimal.json` the RS(6,2) counterexample · `review_gap.py` an independent check of it · `sweep_f16.py`,
`compile_f16.py`, `res_f16.jsonl`, `schemes_f16.json`, `table_f16.json` the search on the GF(16)-based points · `parse_thesis.py` recomputes every scheme printed in the 2024
thesis (`thesis_appendix_recomputed.json` is its output; the thesis itself is not included).
