# A counterexample to the subfield-point comparison for Reed–Solomon repair

**Result.** For RS(6,2) over GF(256), the optimal worst-node bandwidth is **10 bits** on the ISA-L evaluation points and **11 bits** on the prescribed GF(16)-based evaluation points of Dinh et al. This gives a negative answer to the fixed-family comparison posed in [Section IV-B of their ISIT 2022 paper](https://arxiv.org/html/2205.11015#S4.SS2). Their Table IV reported achievable values of 12 bits for both sets. This note gives an impossibility proof and matching constructions.

Consider one erased byte. Each surviving node simultaneously sends a fixed number of bits determined by its own byte and the identity of the failed node. The target must be recovered exactly for every codeword. Binary-linear repair is the model of the cited comparison. The proof below also covers arbitrary deterministic helper functions in this one-round, fixed-bit model.

Use the byte representation of GF(256) with modulus 0x11D, where addition is XOR. Set z = 78, a root of z⁴ + z + 1 = 0 of multiplicative order 15. The paper prescribes a root of this polynomial, and its four choices give equivalent codes under the Frobenius automorphism. In the orders shown, the **exact** bandwidth profiles are

| Evaluation points | Optimal bits for each failed point | Worst |
|---|---|---:|
| 0, 1, 2, 3, 4, 5 | 10, 10, 10, 10, 10, 10 | 10 |
| (0, 1, z, z², z³, z⁴) = (0, 1, 78, 153, 68, 79) | 10, 10, 10, 11, 10, 11 | 11 |

A codeword consists of the evaluations of an affine polynomial f over GF(256). Making the generator systematic preserves the code. Invertible per-coordinate rescaling preserves all bandwidths because a node can undo its own scale locally. For the ISA-L (6,2) Cauchy generator this is explicit: given the two data symbols f(0), f(1), its row at i = 2,3,4,5 stores f(0)/i + f(1)/(i+1) = f(i)/(i(i+1)).

## When ten bits suffice

Let F = GF(256), and fix a failed point t and a reference helper h. For the four other helpers define

$$c_i = (i-h)/(i-t).$$

**Lemma.** Ten bits suffice for binary-linear repair at t if and only if these four elements are linearly independent over GF(2).

**Lower bound.** Write B for the total number of bits and ℓ_j for the contribution of helper j. Given f(j), the target f(t) can still take all 256 values, so the other four helpers must send at least eight bits: B − ℓ_j ≥ 8. Summing over j gives 4B ≥ 40. At equality, each helper sends exactly two bits.

**Necessity.** Put U = f(t), W = f(h). These are independent coordinates for the entire file. Helper i stores

$$f(i)=a_i U+b_i W,\qquad a_i=(i-h)/(t-h),\quad b_i=(i-t)/(h-t).$$

Here b_i denotes an interpolation coefficient; in characteristic two, a_i/b_i = c_i. Write T_i for the binary space of trace multipliers used by helper i: a multiplier v gives the bit Tr(v f(i)). The nondegenerate trace pairing represents every binary-linear functional in this way. At total cost ten, fixing any one helper's stored byte leaves the other four messages in bijection with U. In particular all helper maps have rank two, so each T_i has dimension two.

Fix W = 0. The eight bits from the four non-reference helpers form an invertible linear map of U. If W varies inside the kernel of h's message, its contribution to these eight bits must vanish. Otherwise a nonzero change in U could cancel that contribution, giving the same complete transcript for different targets. Thus every trace with multiplier in b_i T_i vanishes on the kernel of h's message. By trace duality b_i T_i is contained in T_h. Both spaces have dimension two, so they agree.

Consequently the eight target traces must span the direct sum

$$F=\bigoplus_{i\ne t,h}c_i T_h.$$

A binary dependence among the c_i gives a dependence among these subspaces by multiplying it by any nonzero element of T_h. This rules out ten bits for **every** binary-linear repair scheme at such a node.

**Sufficiency.** If the four c_i are independent, let C be their binary span. The ratio set of its 15 nonzero elements has at most 15² = 225 members. Choose β among the 255 nonzero field elements outside this ratio set. Then C and βC intersect only in zero and together span F. The reference helper sends Tr(W) and Tr(βW). Each of the four others sends Tr(v f(i)/b_i) for v = 1, β. XOR with the matching reference bit yields Tr(v c_i U). These eight multipliers span F, so they determine U. This proves the lemma.

## The obstructed nodes and matching eleven-bit repairs

For each t take h to be the first point in the displayed order different from t. On the ISA-L points the four c_i have binary rank four at all six targets. On the prescribed subfield points their ranks in order are 4, 4, 4, 3, 4, 3. The two obstructions can be checked directly:

| t | h | Other helpers i, in order | Corresponding c_i | Binary dependence |
|---:|---:|---|---|---|
| 153 = z² | 0 | 1, 78, 68, 79 | 11, 147, 146, 10 | 11 ⊕ 147 ⊕ 146 ⊕ 10 = 0 |
| 79 = z⁴ | 0 | 1, 78, 153, 68 | 69, 78, 11, 147 | 69 ⊕ 78 ⊕ 11 = 0 |

Eleven bits suffice at both nodes. Let h transmit the three traces of W with multipliers 1, 2, 4. For each other helper use the multiplier list below. For each listed v that helper sends Tr(v f(i)/b_i), and the receiver obtains Tr(v c_i U) by cancellation. An empty list means zero bits.

| t | First helper | Second helper | Third helper | Fourth helper |
|---:|---|---|---|---|
| 153 | 1, 2, 4 | 1, 2, 4 | 2, 4 | empty |
| 79 | 1, 2, 4 | 1, 2, 4 | empty | 1, 2 |

In each row the eight products v c_i have binary rank eight. The three reference bits and the eight other bits therefore give exact repair. Together with the lemma and the ten-bit schemes at the other targets, this proves the profiles.

## Allowing nonlinear helper messages

Suppose ten bits suffice in the deterministic model stated above. The counting argument again forces two bits from each helper. Write φ_i for a helper's encoding. For every fixed W = w, the eight-bit tuple of the four non-reference messages is a bijective function of U. If φ_h(w) = φ_h(w′), those two bijections have the same inverse, the decoder with that reference message fixed. They therefore agree for every U. Because a_i U ranges over F, each φ_i for i ≠ t,h is invariant under translation by b_i(w − w′).

Let H be the additive span of all such differences within the fibers of φ_h. A largest fiber has at least 64 elements and lies inside a coset of H, so |H| ≥ 64. Each other φ_i is constant on cosets of b_i H and takes all four possible values by the bijection above, so |H| ≤ 64. Equality follows. Each φ_i therefore labels exactly the four cosets of b_i H. The reference helper likewise labels exactly the four H-cosets, since each fiber lies in one coset and four labels cover F. Up to relabeling the two-bit outputs, all five helpers send linear quotient maps. The linear obstruction thus also excludes a ten-bit nonlinear repair in this model. For related work on prime-field linearization at cut-set equality under a standing k ≥ 3 assumption, see [Li, Sections II and VII](https://arxiv.org/html/2610.34377v2#S7).

## Reproduce the finite arithmetic

Run `python3 -B check_pointset_gap.py` from this directory. The [verifier](check_pointset_gap.py) computes the multipliers and their ranks and constructs a decoder whose only input is the transmitted bits. It checks all 65,536 affine polynomials at each of six targets on each of the two point sets (786,432 repairs).

The definition of the subfield point set matters. Applying the proved criterion to all 8,008 six-element subsets of GF(16) finds 880 with ten-bit repair at every node. Thus optimizing the choice of *any* subfield six-set also gives ten bits. The strict 10 versus 11 separation concerns the two point-selection rules specified in the 2022 paper.
