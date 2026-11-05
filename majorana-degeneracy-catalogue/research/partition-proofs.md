# Anticoherence at a fixed full Majorana partition

For a partition λ=(λ₁,…,λ_d) of n, in decreasing order, let t_max(n,λ) be the largest anticoherence order attained by **some** pure symmetric state with exactly that complete multiplicity partition. We adjoin order zero as a convention: a maximum of zero means that no first-order anticoherent representative exists. It does not mean that the partition is unrealizable.

A partition generally contains many states, rotation orbits, and SLOCC classes. A verified representative establishes a lower bound only for its own full partition. A largest-point-multiplicity bound establishes an upper bound for every partition containing such a point; it does not establish existence for all those partitions.

## Complete first-order existence criterion

Write g=λ₁. A first-order anticoherent representative exists if and only if

\[
2g<n\qquad\text{or}\qquad
n\text{ is even and }\lambda=(n/2,n/2).
\]

In the strict case, a stronger statement holds: **every** symmetric state in that partition has a first-order anticoherent state in its own orbit under a finite invertible collective filter A^{⊗n}. This preserves all multiplicities and all distinctions between Majorana points. An explicit filtered coefficient vector is not required for the existence conclusion, but should not be claimed to have been tabulated when only this theorem is used.

### Necessity, including the half-size boundary

Rotate a point of multiplicity g to the pole where c_0=⋯=c_(g−1)=0. For the normalized Dicke-index distribution K with probabilities |c_k|², first-order anticoherence requires E[K]=n/2, while K≥g.

If g>n/2 the mean is impossible. If g=n/2, the mean condition forces all probability onto K=n/2. The state is the balanced Dicke state, whose full partition is exactly (n/2,n/2). Thus the other partitions at the half-size boundary are excluded. Conversely, the balanced Dicke state is first-order anticoherent.

This exclusion and the balanced exception are proved in [Baguette, Bastin and Martin, PRA 90, 032314 (2014), PDF p.4](../references/PRA90_032314(2014).pdf#page=4). The source explicitly leaves its remaining classes for further consideration; it must **not** be cited as proving the unrestricted converse below. The 2015 point-group paper proves a sufficiency statement with an additional cyclic-symmetry assumption, rather than the full statement used here: [PRA 92, 052333 (2015), PDF p.7](../references/PRA92_052333(2015).pdf#page=7).

### Sufficiency by an elementary finite-filter argument

Let ψ be any normalized symmetric state for which every Majorana multiplicity is strictly smaller than n/2. Use the spin-axis convention J_z|D_n^(k)⟩=(n/2−k)|D_n^(k)⟩. For V∈SU(2), write V^{⊗n}ψ=∑_k d_k(V)|D_n^(k)⟩, and define

\[
p_+(V)=\sum_{k<n/2}|d_k(V)|^2.
\]

This is positive for every V. Otherwise all coefficients with k<n/2 would vanish. The Majorana polynomial would then have a zero at the pole of multiplicity at least ⌈n/2⌉, contradicting the strict multiplicity hypothesis. Rotations preserve multiplicities.

The function p_+ is continuous on compact SU(2). Since it never vanishes, its minimum ε is strictly positive. This uniform bound is the key step; pointwise positivity alone would not suffice.

Every A∈SL(2,C) has a singular-value decomposition

\[
A=U\operatorname{diag}(e^{s/2},e^{-s/2})V,
\qquad U,V\in\mathrm{SU}(2),\quad s\ge0.
\]

The left unitary does not change the norm. If m_*>0 denotes the smallest positive spin eigenvalue (1/2 for odd n and 1 for even n), then

\[
\|A^{\otimes n}\psi\|^2
=\sum_k |d_k(V)|^2e^{2s(n/2-k)}
\ge\epsilon e^{2s m_*}.
\]

Thus the norm diverges uniformly as s tends to infinity. Its sublevel sets are compact, because s is bounded there and both unitary factors range over compact sets. A finite global minimizer A_0 therefore exists. Its output is nonzero because A_0 is invertible.

Vary A_0 on the left as A(u)=exp(uσ_a/2)A_0, for a=x,y,z and real u. These matrices remain in SL(2,C). For φ=A_0^{⊗n}ψ, differentiation at the minimum gives

\[
0=\left.\frac{d}{du}\|A(u)^{\otimes n}\psi\|^2\right|_{u=0}
=2\langle\phi|J_a|\phi\rangle.
\]

After normalization, φ has all three spin expectations zero, equivalently ρ_1=I_2/2. It is first-order anticoherent. The invertible map A_0 induces a bijective Möbius transformation of the Majorana projective line, so no points merge or split and the complete partition remains λ. Rescaling A_0 to a valid local Kraus operator gives a physical filter with nonzero success probability; it does not change the normalized output.

This is an independently supplied proof. It is consistent with the usual normal-form method: [Verstraete, Dehaene and De Moor, PRA 68, 012103 (2003), Theorem 1 and discussion, PDF p.2](https://arxiv.org/pdf/quant-ph/0105090#page=2), relates local maximally mixed reductions to minimum norm and explicitly discusses infinite filtering limits. The strict root-multiplicity condition is the standard notion of stability for binary forms: [Stoll and Cremona, *On the reduction theory of binary forms*, Definition 4.1, PDF p.7](https://mathe2.uni-bayreuth.de/stoll/papers/reduction-part-1.pdf#page=7). Neither external statement is needed in place of the coercivity proof above.

### Why an orbit closure does not establish the boundary case

At g=n/2, rotate the repeated point so that only coefficients k≥n/2 remain, with c_(n/2)≠0. Applying diag(e^{s/2},e^{−s/2})^{⊗n} and sending s→∞ leaves a nonzero multiple of the balanced Dicke state. If the original partition was not (n/2,n/2), this limit merges other points and changes the full partition. No finite invertible filter can do that. The limit therefore supplies no first-order witness for the original partition. The separate mean argument already excludes it.

## A full-partition restriction from second-order saturation

Suppose a second-order anticoherent state contains a point of multiplicity g. In a pole frame K∈{g,…,n}, and

\[
\mathbb E[(K-g)(n-K)]=\frac{n(n-1-3g)}6\ge0.
\]

When n=3g+1, this expectation is zero. Since the integrand is nonnegative on the allowed support, the Dicke distribution must be supported on {g,n}. The mean n/2 fixes

\[
|c_g|^2=\frac{n}{2(n-g)},\qquad
|c_n|^2=\frac{n-2g}{2(n-g)}.
\]

Both are nonzero. Up to a nonzero scale and phases, the Majorana polynomial is z^g(A+Bz^(n−g)) with AB≠0. Its g-fold pole is the only repeated root: all n−g roots of A+Bz^(n−g) are distinct. Therefore the full partition must be

\[
\lambda=(g,1^{n-g}).
\]

Consequently, at n=7 the other partitions with largest part 2 cannot be second-order anticoherent, and at n=10 the other partitions with largest part 3 cannot be second-order anticoherent. All these partitions are stable, so their maximum order is exactly one. This deduction resolves nine cells beyond the largest-part bound alone. The moment characterization is in [PRA 92, 052333 (2015), PDF p.4, Eqs. (16)–(17)](../references/PRA92_052333(2015).pdf#page=4); the saturation conclusion is independently derived here.

## How the 270 bounds are assembled

1. Enumerate each decreasing integer partition of n exactly once, for n=2,…,12.
2. Apply the complete first-order criterion. Excluded partitions have maximum zero. Every other partition has a proved lower bound of one, whether or not a coefficient vector has been tabulated.
3. Match verified states by their **entire** partition and improve the lower bound to their achieved order.
4. Apply the proved global maximum AC order at that n and the exact moment-polynomial exclusions for g=max λ.
5. Apply the second-order saturation restriction above.
6. Declare an exact maximum only when the resulting lower and upper bounds agree. Otherwise retain the interval without interpreting it as a conjectured optimum.

The global order bounds used are 1,1,2,1,3,2,3,2,3,3,5 for n=2,…,12, backed by the companion exact AC certificates. The input includes 48 verified states; the new tetrahedral representatives add the partitions (3,3,3,3) and (2,2,2,2,1,1,1,1) at n=12. The first has exact maximum two; the second currently has the rigorous range two through three.

The resulting classification has **232 exact maxima and 38 ranges**. Exact counts are 112 at zero, 102 at one, 10 at two, 7 at three, and 1 at five. The ranges are 30 of 1–2, 7 of 1–3, and 1 of 2–3. Every partition through n=7 is classified exactly. No failed search is used as an exclusion.

`verify_partition_bounds.py` independently recomputes the combination from the packaged ingredients, verifies each moment-exclusion polynomial with rational arithmetic, and compares all 270 bounds with `partition-audit.json` and the displayed catalogue. The individual state, multiplicity, and global nonexistence certificates retain their separate exact verifiers.

The subsequent exact constructions in `ac-order-sparse.json` and `ac-order-constructions.json` raise nine lower bounds to two. Six match the existing upper bound, settling their maxima; the three affected eleven-qubit partitions remain between two and three. The independent audit `ac-order-improvements-audit.py` checks exact reduced matrices or moment identities modulo isolated parameter polynomials, and full root multiplicities.
