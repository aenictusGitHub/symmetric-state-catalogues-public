# Anticoherent subspaces

A static, offline catalogue of the largest dimension of a pure anticoherent subspace in the symmetric n-qubit Hilbert space. Open `index.html`; no build or network connection is required. KaTeX, source records and the Majorana viewer are bundled locally.

The grid covers n = 1–24, with selected larger constructions at n = 35, 42, and orders t = 1–7. It is indexed by qubit number n and anticoherence order t. A single value is a proved maximum; an interval gives rigorous lower and upper bounds. Numerical dimensions are shown separately and never silently promoted to rigorous lower bounds. Unresolved grid cells show `Conj. k` and `Num. k` badges when available. If their values agree, only `Conj. k` is shown; proved maxima and proved exclusions omit both badges. The detail panel retains the paper’s Figure 2 values and numerical reports, with confirmed, attained-but-open, or unresolved status. All 33 applicable paper entries, including zero dimensions, remain included. Filters select qubit number, evidence, paper conjectures or cells with numerical reports.

Select a cell to inspect a representative subspace, its exact orthonormal Dicke basis or numerical basis, the orthogonality and cross-reduction checks, exact reconstruction weights and phases when supplied, and the supporting proof/source records. Floating-point checks of decimal components are labelled separately from exact certificates. Numerical source vectors are normalized individually and are only approximately orthogonal to the reported accuracy. Every vector in an AC subspace must be anticoherent: checking the individual basis vectors alone is insufficient. The complete condition is `R_t(|u_a><u_b|) = δ_ab I/(t+1)`.

Each attached basis vector can be viewed on an interactive Majorana sphere with rotation, zoom, symmetry axes and point groups. These constellations and groups describe the selected basis vectors, **not the subspace itself**. A change of basis can change the constellations. Coordinates and geometric group classification are numerical visualizations, separate from the AC certificates.

`data/catalogue.json` is the downloadable dataset; `data/catalogue.js` contains the same records for direct local-file use. The selected-record download includes its representative basis and proof records. Links such as `#n=7&t=2&subspace=...` preserve the selected parameter pair and representative.

Scientific provenance and reproducible calculations are linked from the catalogue's sources, proof records and accompanying research files. The site is ready for static hosting, including GitHub Pages.

## Thirteen qubits: maximum pure-state order is three

`research/n13-t4/index.html` proves `K(13,t)=0` for every `t>=4`, closing the formerly unresolved fourth- and fifth-order cells. A normalized fourth-order AC state would give a trace-one invariant four-copy moment on `Sym^4(V_(13/2))`, of dimension 2380. Its low-multipole contraction equations give 320 independent integer equations on 331 real invariant coordinates. The exact eleven-dimensional kernel is annihilated by the physical trace functional, so every solution has trace zero: a contradiction.

The standalone `research/n13-t4/verify_independently.py` uses only Python's standard library. It reconstructs the highest-weight representations, selected contraction equations, modular ranks, integer kernel and physical trace from the bundled certificate. It covers arbitrary complex states. Neighboring exact N=12 and N=14 witnesses survive the same contraction construction as independent numerical positive controls.

The accompanying cross-catalogue audit checks 125 shared parameter pairs and all 100 local Liège Dicke files. The pure-state catalogue now includes the exact N=14 and N=17 fourth-order states and seven formerly unknown exclusions, including N=13 at orders four and five. Numerical library checks are kept distinct from exact certificates.

## Patterns across all exact maxima

`research/exact-maxima-patterns/index.html` compares every catalogue entry, keeping proved maxima, proved zeros and unresolved intervals distinct. It proves two new universal upper bounds: `7K(N,3) <= N+7` for `N>=6`, and `10K(N,4) <= N+10` for `N>=9`. Their auxiliary dimensions continue the sequence `3,5,7,10`, one larger than the rectangular-tensor counts at orders one through four. Extension to higher orders is unproved. The first improvements over the preceding general bounds occur at N=55 (order three, upper nine becomes eight) and N=99 (order four, upper eleven becomes ten).

The same analysis proves `K(16,6)=0`, closing the previous interval 0–1. Equality in the two-qubit rectangular-tensor bound forces the rank-seven multipole to vanish because its exact channel multiplier is `5 sqrt(69)/8398`. Seventh-order AC would then require twenty orthogonal images in dimension sixteen. Separate physical-tensor and Clebsch–Gordan verifiers confirm the coefficient.

There are also exact links between optimal representatives. The seven-qubit AC2 code projector obeys `P6 = I7 - 2 Tr_1(P7)`, where P6 is the optimal six-qubit AC1 projector. Thus the reduced code state's kernel is the optimal AC1 subspace. Time-reversed vectorization of the centered projector produces the twelve-qubit icosahedral AC5 state. A parallel four-to-eight-qubit construction gives the cubic AC3 state. These require a special maximal-tensor-rank condition; the condition is impossible at source size n>=7, so these links do not establish an infinite family.

All proof sources, exact verifiers, a machine-readable pattern report and the complete dimension table are bundled in that research folder. The general nested-parent constraint `d+s(s+2)k <= N+1` explains why maximum-dimensional AC1 parents need not contain optimal AC2 codes.

The new bounds also make a complete classification of unrounded tensor-packing equality possible through order four. Equality `2k=N` at t=1 occurs only at `(N,k)=(2,1),(4,2),(6,3)`; equality `4k=N+1` at t=2 only at `(7,2)`; equality `6k=N` at t=3 only at `(6,1)`; and equality `9k=N+1` at t=4 never occurs. Thus `6K(N,3)<=N-1` for N>=7 and `9K(N,4)<=N` for N>=4. The proof reduces the unbounded ranges to finite exact channel checks and handles the remaining N=17 case by the existing vectorization obstruction. This concerns unrounded packing equality, not attainment of rounded or stronger bounds.


## General second-order constructions and a new optimal example

`research/t2-general/index.html` supplies explicit analytical prescriptions for arbitrary qubit numbers. Separated triples of Dicke levels give
`K(N,2) >= floor((N+3)/9)` for every `N>=6`. The subsequent prescribed-moment refinement in `research/moment-partitions/index.html` gives
`K(N,2) >= floor((N-1)/8)` for every `N>=17`, with rational weights and recursively constructed radical phases. Combining the two terms yields `max(floor((N+3)/9),floor((N-1)/8))` throughout `N>=6`, since the triple term dominates below seventeen. Exact polynomial certificates cover all eight residue classes, with a separate independent reconstruction. These are constructive lower bounds; a construction attaining the unrestricted maximum for every N remains unproved.

The adjacent-triangle prescription additionally gives an explicit three-dimensional radical basis at N=22. A cyclic-residue search at N=24 produces a five-dimensional basis defined by 23 integer quadratic equations and a unique root in a rational box of radius `10^-60`. Exact contraction proves existence, and independent Dicke splitting verifies all 225 ordered cross entries. The upper bound gives five, proving `K(24,2)=5`.

| N | Previous rigorous interval at t=2 | New interval |
|---:|---:|---:|
| 22 | 2–5 | 3–5 |
| 24 | 3–5 | 5 |
| 35 | 4–8 | 5–8 |
| 42 | 5–9 | 6–9 |

All earlier representatives remain available. Each new representative includes analytical coefficients and certified Majorana coordinates. The general cyclic ansatz is a search prescription, with no general feasibility or optimality guarantee.

```bash
python3 research/t2-general/separated_triples.py
python3 research/t2-general/verify_separated_independently.py
python3 research/t2-general/four_bands.py
python3 research/t2-general/fourband-independent-audit.py
python3 research/t2-general/construct_exact.py --n 40 --exact-json /tmp/ac40-exact.json
python3 research/t2-general/adjacent-triples.py
python3 research/t2-general/n24-k5/verify.py
python3 research/t2-general/n24-k5/verify_independently.py
```

## Prescribed-moment convex partitions and improved radical constructions

`research/moment-partitions/index.html` formulates AC2 constructions as convex combinations of the lifted Dicke levels `(m,m²)` at the prescribed target `(0,N(N+2)/12)`. Disjoint blocks on globally three-separated levels give exact square-root coefficients. An arbitrary common Tverberg intersection is insufficient: its moments must equal the AC target. A Fan–Pall graph construction allows overlapping basis supports and provides an exact radical prescription for the corresponding joint compression.

The new four-band moment weights preserve those prescribed moments and admit explicit quadrilateral phase cancellation. They prove `K(N,2)>=floor((N-1)/8)` for every `N>=17`, improving the previous four-band term by one on seven of the eight residue classes. The proof uses exact polynomial inequalities for all sizes, independently audited; numerical checks only verify the rendering. For example the combined general construction increases from 4 to 5 at N=41, from 5 to 6 at N=49, from 6 to 7 at N=57, and from 7 to 8 at N=65. These compare general construction formulas, not best size-specific published dimensions.

The catalogue also includes simpler real radical representatives at `(N,k,t)=(12,2,2)` and `(22,3,2)`. All 117 ordered two-qubit cross entries are verified exactly. Their Majorana representations are available as for the other bases. Existing grid intervals are unchanged: stronger size-specific witnesses already cover every displayed cell to which the new family applies.

A separate exact limitation proves `k<=ceil(N/9)` for any AC2 subspace supported on globally three-separated Dicke levels, including overlapping vectors and nonuniform supports. This is a bound on the sparse construction class, not the unrestricted maximum. The new four-band construction admits nearby levels and cancels their cross terms.

```bash
python3 research/moment-partitions/families_fourband.py
python3 research/moment-partitions/search_independent_family_audit.py
python3 research/moment-partitions/search_exact_examples.py
python3 research/moment-partitions/geometry_verify.py
```

## QEC and joint-compression upper-bound investigation

`research/t2-qec-upper/index.html` records a further investigation using intrinsic quantum-code weight enumerators, shadow positivity, and unrestricted Hermitian ancillary coefficients. No stronger unrestricted upper bound was obtained. Exact rational points independently verified via factorial six-j sums show that the specified expanded enumerator relaxation admits the current candidate dimensions `(N,k)=(35,8),(41,9),(42,9)`. These points are not physical AC subspaces and do not establish existence. Finite matrix-witness searches produced no exclusion certificate.

A proved structural consequence is `5k<=N+5-g` when all binary forms in an AC2 subspace share a factor of degree g. It follows from injectivity of multiplication by homogeneous degree-four binary forms. This is conditional: g can be zero. In particular, a code contained in one Dicke-parity sector at odd N has `5k<=N+4`; this does not cover arbitrary codes or codes using both sectors. No catalogue values changed.

```bash
python3 research/t2-qec-upper/verify_lp_feasible.py
```

## Upper bounds from nesting: what is proved

`research/t2-nesting-upper/index.html` gives a nesting derivation of the existing unrestricted upper bound `min(floor((N-2)/4),floor((N+5)/5))` for `N>=8`. Normalized creation and annihilation images of any AC2 code are structured `2k`-dimensional spaces on `N+1` and `N-1` qubits. Their prescribed spinor compressions, combined with eigenvalue interlacing, recover the two inequalities. These image spaces are not AC1 subspaces. The derivation is universal; 418 exact physical operator checks and an exact N=18 example independently check its implementation.

For a same-N AC1 parent of dimension d containing an AC2 code of dimension k, the exact dimension count is `d+3k<=N+1`. The first-order maximum is an upper bound on d and cannot be substituted without an extension theorem. The known N=18 example has largest containing parent dimension five, although `K(18,1)=7`. No stronger unrestricted upper bound was established by this investigation; no catalogue interval changes.

```bash
python3 research/t2-nesting-upper/verify_conditional_parents.py
```

## Stronger lower bounds from varying the first-order parents

`research/t2-nesting-bounds/index.html` proves the universal radical construction `K(N,2)>=floor(N/8)-1` for every `N>=40`, improving the earlier four-band term by two dimensions. Exact polynomial coefficient certificates cover every qubit number; an independent audit reconstructs the moment weights and all phase-closure inequalities. Each constructed k-dimensional AC2 space lies inside an explicit 2k-dimensional AC1 parent. At that stage, combining with the separated triples gave `max(floor((N+3)/9),floor(N/8)-1)` for `N>=6`. The later prescribed-moment construction above strengthens the second term to `floor((N-1)/8)`.

Variable first-order parents in a fixed Dicke-parity sector also give exact algebraic witnesses of dimension five at N=35 and six at N=42. The catalogue intervals improve from 4–8 to 5–8 and from 5–9 to 6–9. Rational contraction certificates isolate roots of 125 and 180 quadratic equations in radius `10^-50` boxes; a separate verifier checks all 225 and 324 ordered two-qubit cross-reduction entries. Coefficients have exact implicit algebraic definitions, with decimal previews for visualization. The Majorana multiplicities are certified with rational root disks. These lower bounds do not establish maximality, and no unrestricted upper bound changes.

```bash
python3 research/t2-nesting-bounds/lower_fourband_cert_all.py
python3 research/t2-nesting-bounds/lower_fourband_independent_audit.py
for example in n35-k5 n42-k6; do
  python3 research/t2-nesting-bounds/verify_parity_certificate.py "research/t2-nesting-bounds/$example/certificate.json"
  python3 research/t2-nesting-bounds/verify_parity_independently.py "research/t2-nesting-bounds/$example"
done
```

## Constructing second-order subspaces inside first-order parents

`research/t2-nested-t1/index.html` gives the exact reduction to five compressed quadrupole matrices inside any analytical first-order parent. A separated reflection-module subfamily reduces further to one diagonal moment condition, solved by explicit low/high pairings with radical weights. The displayed N=18 example has dimension two; all 36 two-qubit cross-reduction entries are verified exactly. An inertia certificate proves that two is the maximum within the specified seven-dimensional reflection-pair parent. The corresponding nine-dimensional parent at N=24 allows dimension at most three.

The same note proves an exact extension criterion for any second-order code C of dimension k: added first-order vectors must form a first-order subspace in the complement of C, JxC, JyC and JzC. For the catalogue's optimal four-dimensional N=18 code, that complement has dimension three and allows exactly one added dimension. Its largest first-order parent therefore has dimension five, even though K(18,1)=7. Two independent exact verifiers check the projector, spin compressions, and an explicit five-dimensional extension. Searching only maximum-dimensional first-order parents would miss this code.

That earlier note gives construction and parent-specific results. The subsequent variable-parent improvements are described above.

```bash
python3 research/t2-nested-t1/nested_examples.py
python3 research/t2-nested-t1/extension_obstruction.py
python3 research/t2-nested-t1/extension_independent_audit.py
python3 research/t2-nested-t1/verify_reflection_inertia.py
```

## Why the cyclic prescription cannot saturate infinitely often

`research/t2-saturation/index.html` proves that for every `N>=420`, no second-order AC subspace of dimension `k=floor(N/5)+1` can be invariant under a collective rotation by `2*pi/k` about any axis. On the subsequence `N=5k-1`, exclusion already holds for `k>=62`.

The proof first shows that every cyclic character sector must have dimension at most one, using an endpoint-variance argument. The error-space projector inequality then supplies a seven-residue weighted moment certificate and an exact contradiction. Two separate verifiers check the identities and threshold gaps. Thus the obstruction applies to every cyclic-invariant candidate of the target dimension, not only a particular choice of coefficients or a smooth ansatz.

This is a restriction on symmetry, **not an improved unrestricted bound** on `K(N,2)`. All catalogue dimensions and intervals remain unchanged. No infinite family attaining the unrestricted general upper bound has been established. The research note keeps numerical search failures separate from the proof.

```bash
python3 research/t2-saturation/construction_cyclic_exclusion.py
python3 research/t2-saturation/cyclic-exclusion-independent-audit.py
```

## A stronger general second-order upper bound

`research/t2-spin2-inertia/index.html` proves `5k≤N+5` for every `N≥4` and AC order `t≥2`. With auxiliary spin 2 and `X=J·S`, the operator `W=X²+(j+3)X−2j(j+1)I` has zero compression on the code tensored with the five-dimensional auxiliary space. Its only positive sector is total spin `j+2`, of dimension `N+5`; projection into that sector proves the bound.

Combining this with the spin-3/2 bound gives `K(N,2)≤min(floor((N−2)/4),floor((N+5)/5))` for every `N≥8`. The first integer improvement is at N=34, where the upper bound becomes seven instead of eight. Thus `floor((N−2)/4)` cannot be the unrestricted maximum for all sizes. A general construction attaining the combined bound has not been established.

When introduced, this bound improved the N=42, t=2 interval from 4–10 to 4–9; the later general construction raises its lower bound to five. Separately, `research/n20-t2/index.html` proves `K(20,2)=4` using a four-dimensional algebraic construction. The later exact reconstructions at N=18,19 prove dimension four there as well. Every earlier representative is retained. Exact constructions now prove the second-order maxima throughout 8≤N≤20, and at N=24.

```bash
python3 research/t2-spin2-inertia/verify.py
python3 research/t2-spin2-inertia/verify_independently.py
```

## Twenty qubits: a certified four-dimensional second-order subspace

`research/n20-t2/basis.json` explicitly defines an orthonormal basis using the unique root of 72 rational quadratic equations in a rational box of radius `10^-35`. Twelve further coordinates are fixed rationally. Both verifiers reconstruct the complete first- and second-moment equations, including every cross term. Exact integer contraction inequalities prove existence. These are algebraic coefficients specified by an isolating system, rather than an asserted radical formula.

The earlier spin-3/2 bound gives `k≤4` at N=20, so the certified construction proves `K(20,2)=4`. Each basis vector has separately certified Majorana root disks and a numerical sphere rendering. The original numerical example and the earlier two-dimensional Molien witness remain available.

```bash
python3 research/n20-t2/verify.py
python3 research/n20-t2/verify_independently.py
python3 research/n20-t2/verify_geometry.py
```

## First-order optimality for every qubit number

`research/t1-optimality/index.html` proves the exact formula `K(N,1)=floor(N/3)+1` for every `N≥4`. Coupling the physical spin `j=N/2` to an auxiliary spin 1 gives `W=J·S`, whose eigenvalues are `j, −1, −(j+1)` with multiplicities `N+3, N+1, N−1`. Every first-order AC subspace produces zero compression on its tensor product with the three-dimensional auxiliary space. Projection into the positive eigenspace is injective, so `3k≤N+3`. The existing analytical bases attain the resulting bound for every size.

The eleven remaining first-order intervals are now exact maxima: N=14,15,16,17,18,19,20,22,24,35,42. No higher-order bound changes. All paper conjectures represented in the catalogue are now confirmed. The small cases remain `K(1,1)=0` and `K(2,1)=K(3,1)=1`. The sequence agrees with OEIS A008620 from N=4 onward under the convention `a(N)=floor(N/3)+1`.

The proof is general and applies to arbitrary complex subspaces. Symbolic identities and independent physical-matrix checks validate its implementation; no numerical infeasibility claim is used.

```bash
python3 research/t1-optimality/verify.py
python3 research/t1-optimality/verify_independently.py
```

## Molien constructions and six stronger lower bounds

`research/molien/index.html` explains the invariant and single-character polynomial constructions for quaternion, tetrahedral, octahedral and icosahedral symmetry. All 65 nonzero constructions in seven families are supplied in `research/molien/bases.json`, with exact polynomial spanning bases, orthogonalized polynomials, normalized analytical Dicke coefficients, and achieved AC orders. The six representatives that strengthen existing lower bounds are also selectable in the catalogue with Majorana spheres.

| N | t | Previous interval | Updated interval |
|---:|---:|---:|---:|
| 20 | 2 | 1–4 | 2–4 |
| 22 | 3 | 0–3 | 1–3 |
| 24 | 2 | 1–5 | 3–5 |
| 24 | 3 | 1–3 | 2–3 |
| 42 | 2 | 1–10 | 4–10 |
| 42 | 3 | 1–6 | 2–6 |

The octahedral sign character supplies the 22-qubit third-order example. Each character sector is treated separately; direct sums of different characters are not asserted to be AC. Higher achieved orders are checked exactly, and inherited lower orders are included. At the Molien step no upper bound changed and none of the six maxima was settled. The later algebraic construction in `research/n20-t2/` now proves `K(20,2)=4`.

The polynomial construction verifies every reduced cross matrix with exact arithmetic. An independent verifier instead reconstructs compressed spin moments and character dimensions. The normalized Dicke convention includes the Majorana sign and binomial factor explicitly. Molien dimensions count these constructions rather than unrestricted maxima.

```bash
python3 research/molien/construct.py
python3 research/molien/verify_independently.py
python3 scripts/build_data.py
python3 research/molien/build_page.py
python3 research/verify_catalogue_independently.py
```

`research/molien/improvements.json` records every changed entry against the saved pre-Molien snapshot. Existing historical improvement reports retain the bounds at the time those results were added.

## Exact results added from ACLiege

Three sparse numerical support patterns now have independently verified radical constructions: `(N,k,t) = (12,2,2), (13,2,3), (18,2,3)`. The middle case attains the complementary Schmidt-space bound and proves `k_max(13,3)=2`. The twelve-qubit case is now also settled: the rank-one complement certificate proves `k_max(12,2)=2`. The general third-order saturation obstruction now proves `k_max(18,3)=2`. These constructions preserve the source support patterns; they are not asserted to equal the original decimal subspaces.

The tensor-Gram consequence `k(floor(t/2)+1)^2 <= N+1`, combined with the paper's Eq. E2, also gives `k_max(10,2)=2`. No claim of literature priority is made for this bound.

Different records may describe equivalent subspaces or different bases of the same subspace. The independent audit checks every exact cross-reduction matrix entry, along with normalization, orthogonality, dimension bounds, paper conjectures and evidence status. Current counts are recorded in `research/catalogue-independent-audit.json`.

## Reproduce

Python 3, NumPy and SymPy are required only to regenerate or verify the data:

```bash
python3 scripts/build_data.py
python3 research/verify_catalogue_independently.py
python3 research/acliege/radical-subspaces/verify.py
```

The builder uses the bundled checked source snapshots and a geometry cache. Remove `research/geometry-cache.json` to recompute the cached Majorana coordinates. The N=16 analytical basis uses a separate exact square-free certificate, avoiding expensive coefficient-field construction; verify it with `python3 research/analytical-constructions/verify_n16_geometry.py`, or regenerate its coordinates with `--write`. Numerical geometry never serves as an anticoherence certificate. The supplied Julia file is preserved under `references/`; it is never executed by these scripts.

Detailed source checks are in `research/paper/verification.json`, `research/paper/radical-independent-verification.json`, and `research/catalogue-independent-audit.json`. The source notes identify malformed Julia vectors and caveats in the printed dimension formulas, rather than importing them as valid witnesses.

## General bounds for all qubit numbers and AC orders

`research/general-bounds/index.html` proves the rectangular-tensor family `k m_h(m_h+h) ≤ N−h+1`, where `m_h=floor((t−h)/2)+1` and `0≤h≤t`. This contains the Schmidt and ordinary tensor-Gram bounds. Optimizing gives `k(r+1)(r+2)≤N` for odd order `t=2r+1`; for even order `t=2r≥2`, it gives `k(r+1)^2≤N+1` together with the necessary nonzero-existence condition `N≥(r+1)^2`. In particular, pure states require `(t+2)^2≤4N+1`.

General exact channel tests exclude equality unless suitable integer spin ranks have zero multiplier. At third order, saturation `6k=N≥12` requires `2L(L+1)=(N−1)(N+2)`. The odd-order proof uses a general projector-vectorization obstruction, with no size-specific certificate. A further pure-state equality test excludes `N=20,t=7`.

All 147 catalogue entries are evaluated with these bounds, while stronger existing certificates are retained. Thirty upper bounds improve; thirteen entries become settled. Newly proved positive maxima are `k_max(10,3)=k_max(11,3)=k_max(12,3)=1`, `k_max(18,3)=2`, and `k_max(18,5)=k_max(20,5)=k_max(24,5)=1`. Six further entries become proved zero exclusions. Full changes and exact equality metadata are in `research/general-bounds/improvements.json`.

```bash
python3 research/general-bounds/verify_rectangular.py
python3 research/general-bounds/verify_odd_saturation.py
python3 research/general-bounds/verify_t3_saturation.py
python3 research/general-bounds/verify_general.py
```

The arguments apply to arbitrary complex subspaces and every admissible size. Finite exact checks validate the formulas and implementation; they do not replace the general proof. No claim of literature priority is made.


## Sharper second-order bound from an auxiliary spin

`research/ancilla-inertia/index.html` proves `4k≤N−2` for every `N≥8` and every AC order `t≥2`. With auxiliary spin 3/2 and `X=J·S`, the operator `W=X²−(N−4)X/2−5N(N+2)I/16` has zero compression on the code tensored with the four-dimensional auxiliary space. Its spectrum has `N−2` positive directions, `3N+6` negative directions and no zeros. Projection into the positive subspace gives the bound.

Five entries improve: N=13 at order two is now exactly 2; N=16 and N=17 are exactly 3; the upper bounds at N=20 and N=24 fall to 4 and 5. Existing exact witnesses attain this general bound at every N=8 through 17. The boundary N=7 has zero eigenvalues and is explicitly excluded from the theorem's domain, preserving its known dimension-two code.

```bash
python3 research/ancilla-inertia/verify_general.py
python3 research/ancilla-inertia/verify_independently.py
python3 research/n13-t2/inertia/verify.py
```

The analytic proof covers all admissible sizes; independent symbolic and exact physical-matrix checks validate the formulas and conventions. No claim of literature priority is made.

## Twelve-qubit second-order maximum

`research/n12-t2/index.html` gives a complete proof that a three-dimensional second-order AC subspace cannot exist at N=12. The existing exact two-dimensional basis makes the bound sharp. No symmetry, reality, or sparse-support assumption is used.

The contradiction is the exact identity `0 = 9015 + 1031316 a3 + 3936270 p0 + 20020 p2 + 728365 p6`, where every variable is nonnegative. Two independent verifiers reconstruct the recoupling coefficients and check the rational identity:

```bash
python3 research/n12-t2/verify.py
python3 research/n12-t2/verify_independent.py
```

The complete derivation is also downloadable as `research/n12-t2/proof.tex`; no claim of literature priority is made.


## Analytical constructions attaining the paper’s dimensions

`research/analytical-constructions/index.html` gives a general first-order construction of dimension `floor(N/3)+1` for every `N≥4`. Reflection pairs and explicit central blocks attain every positive first-order value in Figure 2, with the established small-size examples at N=2 and N=3. The selected N=7 case now has an exact dimension-three basis and a matching upper bound. The auxiliary spin-1 bound now proves this entire first-order family optimal for every N≥4.

At second order, exact bases give dimension two at N=11 and dimension three at N=16. The auxiliary-spin inertia theorem now proves `k_max(16,2)=3`; the analytical construction is optimal. Earlier spin-error saturation arguments remain valid at N=11,15,19,35. The dimension-three cases at N=14,15,17 also have exact certificates and matching upper bounds; see the construction section below.

Reproduce these exact checks with:

```bash
python3 research/analytical-constructions/verify_t1.py
python3 research/analytical-constructions/verify_t2.py
python3 research/analytical-constructions/verify_gram_saturation.py
python3 research/verify_displayed_formulas.py
```

The constructions are independently derived from moment and cross-term equations; they are not asserted to equal the original decimal subspaces or to be new in the literature.


## Three-dimensional second-order examples at N=14,15,17

`research/algebraic-k3/index.html` contains exact constructions for all three remaining proposed dimension-three cases. N=14 has coefficients in nested square roots. N=15 and N=17 use uniquely isolated real algebraic roots of explicit rational polynomial systems; no simple radical expressions are claimed for those two cases. Their rational contraction certificates prove existence and uniqueness, and every full-span cross reduction follows exactly from the defining equations. The downloaded basis records include the complete parameter definitions and rational isolating intervals.

The witnesses and the auxiliary-spin inertia bound prove `k_max(14,2)=k_max(15,2)=k_max(17,2)=3`. All proposed dimensions in Figure 2 have exact attainment or a proved zero exclusion; cells are labelled optimal only when their lower and upper bounds agree.

```bash
python3 research/algebraic-k3/n14/verify.py
python3 research/algebraic-k3/n15/verify.py
python3 research/algebraic-k3/n17/verify.py
python3 research/algebraic-k3/verify_geometry.py
python3 research/algebraic-k3/verify_certificates_independently.py
```

Majorana coordinates for these bases use rational interval enclosures and disjoint Rouché disks to certify simple finite roots. The rendered coordinates are approximations. Named algebraic-root coefficients must be interpreted together with each record’s `algebraicDefinition`; they are not free parameters or rounded exact amplitudes.


## Sharp eight-qubit dimensions and a general first-order bound

`research/n8/index.html` proves the complete eight-qubit classification: `k_max(8,1)=3`, `k_max(8,2)=k_max(8,3)=1`, and no nonzero subspace at orders four or higher. The exclusions apply to arbitrary complex subspaces.

For first order, saturation of `k=N/2` forces the projector to have only scalar and maximal-rank spin-tensor components. Time-reversed vectorization then produces a pure symmetric `2N`-qubit state of order `N−1`. At N=8, the exact identity `0=100+136 p0+115 a9` excludes the auxiliary sixteen-qubit state. For every larger even N, the tensor-Gram dimension bound already excludes it. Thus `k_max(N,1)≤N/2−1` for every even N≥8. The existing witnesses also make this bound sharp at N=10 and N=12, with dimensions four and five.

For second order at N=8, a hypothetical dimension-two subspace leaves a pure complement to its four spin-error images. Its exact identities force `0=5831+217971 a5`, contradicting nonnegativity of the squared multipole norm.

```bash
python3 research/n8/t1/verify.py
python3 research/n8/t1/verify_independently.py
python3 research/n8/t2/verify.py
```

The verifiers check rational recoupling identities, the maximal-spin vectorization, partial traces, and exact lower witnesses. The separate first-order verifier reconstructs the required tensor contractions directly using 3j coefficients and highest-weight vectors.

## Nine qubits: second-order maximum dimension is one

`research/n9-t2/index.html` proves `k_max(9,2)=1` for arbitrary complex subspaces. A putative orthonormal pair gives a decomposable antisymmetric vector `w=u∧v`. Its whole-subspace cross conditions imply matrix contraction equations for the rotational average of `|w⊗w><w⊗w|`, supported on the 825-dimensional Young-symmetry space of shape `(2,2)`.

The exact certificate has 169 real invariant coordinates, 156 independent rational equations, and a 13-dimensional kernel. A strictly positive integer combination of all 57 diagonal entries of the multiplicity blocks vanishes on that kernel. Positivity therefore forces every block to vanish, contradicting trace one. The existing exact pure state attains dimension one, so the paper's proposed value is now confirmed. Numerical infeasibility is not used as proof.

```bash
python3 research/n9-t2/verify.py
python3 research/n9-t2/verify_independently.py
```

The second verifier reconstructs the physical tensor-word contractions separately. The bundled LaTeX proof specifies the basis, contractions and positive identity; `certificate.json` contains all integer weights. These are catalogue calculations, with no claim of literature priority.

## Thirteen qubits at second order: maximum dimension two

`research/n13-t2/index.html` proves `k_max(13,2)=2` for arbitrary complex subspaces. Coupling spin 13/2 to an auxiliary spin 3/2 gives a Hermitian operator with eleven positive and forty-five negative eigenvalues, no zeros, and zero compression to the `4k`-dimensional code–auxiliary space. Projection into its positive eigenspace is injective, so `4k≤11`. This excludes dimension three without any symmetry or sparsity assumption. The existing exact two-dimensional family attains the bound. Full LaTeX and a standard-library-only exact physical verifier are in `research/n13-t2/inertia/`.

The sparse basis extends to an exact phase family. Distinct phases in `[0,π]` give inequivalent subspaces under collective rotations, as proved by a time-reversal invariant. Every member now attains the proved global maximum. Run `python3 research/n13-t2/verify_family.py` for the symbolic full-span, determinant and invariant checks.

The earlier conditional certificate in `research/n13-t2/containing-u/` excludes codes containing a two-component Dicke state. It remains valid but is superseded by the unrestricted inertia theorem. Its exact certificate and verifier are retained for reproducibility.

The historical `research/n13-t2/projector-moment/` relaxation has a numerically consistent candidate with solver status `optimal_inaccurate`. The global nonexistence theorem confirms that this relaxed moment cannot represent an actual dimension-three AC code. Its generator, optimizer-free candidate check and known-code controls are retained; they are not needed for the new proof.

The exact scalar relaxation admits a rational feasible point. Its portable checker reconstructs all 450 inequalities and nine equality rows; the lower-bound checker verifies all 36 whole-subspace cross entries. Stronger matrix and partial-transpose calculations remain exploratory numerical evidence.

```bash
python3 research/n13-t2/verify_lower.py
python3 research/n13-t2/verify_scalar.py
```

## Eleven qubits: first-order maximum dimension is four

`research/n11-t1/index.html` proves `k_max(11,1)=4`, confirming the paper's proposed dimension. A hypothetical five-dimensional subspace gives a normalized exterior product in the 792-dimensional space `Λ⁵ V_(11/2)`. Its full compressed-spin conditions yield 100 independent rational equations on 111 invariant coordinates. A strictly positive integer identity on all 49 block diagonals forces every positive block to vanish, contradicting normalization. The argument applies to arbitrary complex subspaces and uses no numerical infeasibility claim.

The existing exact four-dimensional Julia basis supplies the lower bound. Its simplified displayed coefficients and all 64 ordered one-qubit cross-reduction entries are checked independently. The complete LaTeX proof, integer certificate, basis, equations and portable verifiers are bundled with the research page.

```bash
python3 research/n11-t1/verify.py
python3 research/n11-t1/verify_independently.py
python3 research/n11-t1/verify_lower.py
```

## Thirteen-qubit first-order maximum

`research/n13-t1/index.html` proves `k_max(13,1)=5`. A hypothetical dimension-six subspace gives a normalized vector in the 3003-dimensional exterior space `Λ⁶ V_(13/2)`. Its full compressed-spin conditions imply 716 independent exact equations on 735 invariant coordinates. A strictly positive integer identity on all 151 block diagonals forces the positive rotational average to vanish, contradicting normalization. Two separate verifiers reconstruct the physical equations and check the certificate for arbitrary complex subspaces. The existing five-dimensional analytical central-block basis attains the bound. This finite-size certificate is retained as an independent proof. The new auxiliary spin-1 argument in `research/t1-optimality/index.html` proves the first-order family `floor(N/3)+1` optimal for every `N≥4`.

```bash
python3 research/n13-t1/verify.py
python3 research/n13-t1/verify_independently.py
python3 research/n13-t1/verify_lower.py
```

## Fourteen qubits at fourth order: radical states and optimal dimension

`research/n14-t4/index.html` gives exact five-component Dicke expressions for the numerical state in `ACLiege.jl`, and for its rotation-inequivalent complex conjugate. Both have proper rotational symmetry group D3. Two independent exact checks give `rho_4=I_5/5` and fifth-reduction purity `15209/90750`, so their AC order is exactly four. The tensor-Gram upper bound `9k≤15` proves `K(14,4)=1`. This settles the former rigorous interval 0–1 while retaining the original numerical representative. The later exact three-copy positivity certificate in `research/n14-t5/index.html` excludes fifth order, proving that four is also the globally maximum pure-state AC order.

```bash
python3 research/n14-t4/verify.py
python3 research/n14-t4/verify_independently.py
```


## Sixteen qubits at third order: an optimal analytical plane

`research/n16-t3/index.html` gives an explicit two-vector basis with radical unit phases. All 64 ordered cross-reduction entries satisfy `R_3(|u_a><u_b|)=delta_ab I_4/4`, so every normalized complex superposition is third-order AC. The rectangular-tensor bound `6k≤16` proves `K(16,3)=2`, closing the former interval 1–2. Each displayed basis vector has exactly AC order three. A common collective rotation and separate basis phases match the supplied numerical Julia basis within `5.10e-10`; the symbolic proof is independent of that approximate match.

```bash
python3 research/n16-t3/verify.py
python3 research/n16-t3/verify_independently.py
python3 research/n16-t3/verify_geometry.py
```

## Seventeen qubits at fourth order: analytical states and optimal dimension

`research/n17-t4/index.html` constructs a six-component radical state with proper rotational symmetry D3 and its rotation-inequivalent complex conjugate. Both have `rho_4=I_5/5` and fifth-reduction purity `5443/32487>1/6`, hence exactly AC order four. These are new representatives using the numerical Julia state's support; equality or rotation equivalence to that source is not asserted. The two branches define separate one-dimensional subspaces, rather than a two-dimensional AC span.

A hypothetical dimension-two fourth-order AC code would saturate `9k=18`. The rank-two tensor channel then forces its projector to contain only scalar and rank-16 multipoles. Its Schmidt-partner projector on sixteen qubits has rank four and the same multipole support. Time-reversed vectorization of its traceless part would yield a pure symmetric 32-qubit state of AC order 15. This requires 64 nonzero orthogonal spin-tensor images in dimension 33, a contradiction. Thus `K(17,4)=1` for arbitrary complex subspaces, closing the former interval 0–2. The complete proof and independent exact recoupling, partial-trace and time-reversal checks are bundled. Existence of a different fifth-order pure state is not settled.

```bash
python3 research/n17-t4/verify.py
python3 research/n17-t4/verify_independently.py
python3 research/n17-t4/verify_geometry.py
```


## Seventeen qubits at third order: a radical family of optimal planes

`research/n17-t3/index.html` gives an explicit two-parameter radical construction with two five-component Dicke vectors. The catalogue displays a simple member with rational probability lists `(95,159,229,125,160)/768` and `(146,181,145,215,81)/768`. Every ordered three-qubit cross reduction is exactly `delta_ab I_4/4`; an independent physical-spin calculation verifies the same entire complex span. The rectangular bound `6k≤17` proves `K(17,3)=2`, closing the former interval 1–2.

A second exact rational-parameter choice in the same family reconstructs the reported Julia vectors within `2.5e-10` after a common collective rotation and separate basis phases. Its full radical expressions are in `source-reconstruction.json`. `verify_family.py` checks the general moment, cross-product and phase identities and certifies the source-near parameter domain using rational intervals. Equality to rounded decimals or a classification of all optimal planes is not asserted.

```bash
python3 research/n17-t3/verify.py
python3 research/n17-t3/verify_independently.py
python3 research/n17-t3/verify_family.py
python3 research/n17-t3/verify_geometry.py
```

## Eighteen qubits at second order: an explicit radical four-dimensional basis

`research/n18-t2/index.html` gives four radical basis vectors: reflection-even and reflection-odd vectors on Dicke residue zero modulo three, and a reflected pair on residues one and two. Exact verification includes all 144 complex cross-reduction entries. Hence every normalized complex superposition is second-order AC. The general auxiliary-spin bound `4k≤N−2` proves `K(18,2)=4`, closing the former interval 2–4. The numerical source is retained separately; a common rotation and a change of basis approximately identify its span with the exact reconstruction.

```bash
python3 research/n18-t2/verify.py
python3 research/n18-t2/verify_independently.py
python3 research/n18-t2/verify_geometry.py
```

## Nineteen qubits at second order: an exact algebraic four-dimensional basis

`research/n19-t2/index.html` specifies four sparse vectors through a unique isolated real root of eighteen integer quadratic equations. This is an exact algebraic basis, rather than a short formula in radicals. All variables have explicit rational isolation intervals. A rational approximate inverse proves that the fixed-point map is a strict contraction and maps the box strictly into itself. The primary verifier reconstructs the complete physical moment system and exact polynomial dependencies; the independent verifier instead builds all 144 ordered Dicke cross reductions directly from the exported coefficient expressions and proves they vanish at the certified root.

The bound `4k≤17` excludes dimension five, proving `K(19,2)=4` and closing the former interval 2–4. Two residue vectors and their unitary coefficient reversals make up the basis; no complex conjugation is used in this reversal. The approximately corresponding source plane remains a separate numerical record. Exact rational Rouché disks certify the Majorana multiplicities of all four representatives.

```bash
python3 research/n19-t2/verify.py
python3 research/n19-t2/verify_independently.py
python3 research/n19-t2/verify_geometry.py
```

## Fourteen qubits: fifth order excluded

`research/n14-t5/index.html` proves `K(14,t)=0` for every `t>=5`. The invariant three-copy moment has 47 real coordinates, and 45 independent exact contraction equations leave two parameters. Three nonnegative multiplicity diagonals have a strictly positive relation; their coefficient forms span both parameters, so positivity forces the entire moment to zero, contradicting trace one. The independent standard-library verifier rebuilds all equations and checks the integer certificate. This applies to arbitrary complex pure states. The known fourth-order radical witness therefore has globally maximum AC order four.

## Fourteen qubits: the complete dimension row is exact

`research/n14-t3/index.html` proves `K(14,3)=1` for arbitrary complex subspaces. A hypothetical plane produces a normalized positive four-copy moment on the 4200-dimensional Young space of shape (2,2). An exact rational combination of 939 contraction equations gives 25 positive-definite multiplicity blocks; all 190 LDL pivots are checked exactly. The independent physical-word verifier uses only the Python 3 standard library. The equations force the resulting strictly positive functional to vanish, a contradiction. The elementary state `(2|D14^0> + sqrt(13)|D14^7> + 2|D14^14>)/sqrt(21)` supplies the lower bound.

Together with the independent fifth-order exclusion and existing exact constructions, the entire N=14 row for orders t=1,…,7 is now `5,3,1,1,0,0,0`.

## Fifteen qubits: third-order dimension one is optimal

`research/n15-t3/index.html` proves `K(15,3)=1` for arbitrary complex subspaces. A normalized four-copy moment of a hypothetical AC3 plane lies in the 5440-dimensional Young-(2,2) space. An exact rational combination of 1197 necessary equations gives 27 strictly positive-definite blocks with 230 positive rational LDL pivots, forcing that moment to vanish. The standard-library verifier reconstructs every equation and checks the elementary four-component lower witness. The same page gives a short saturation proof of the already established `K(15,2)=3`, with links to the exact algebraic basis.

## Fifteen-qubit fourth-order exclusion

The exact result is **(K(15,2), K(15,3), K(15,4)) = (3, 1, 0)**. An independently reconstructed four-copy partial-transpose positivity certificate excludes every pure AC4 state, with arbitrary complex coefficients. All higher orders are excluded by nesting. See [the proof and certificate](research/n15-t4/index.html), [LaTeX proof](research/n15-t4/proof.tex), and [independent exact verification](research/n15-t4/independent-verification.json). Both the pure-state and subspace catalogues are updated.

## Seventeen qubits: fifth and sixth order excluded

The complete N=17 row is now **6, 3, 2, 1, 0, 0, 0** for orders one through seven. The [fifth-order proof](research/n17-t5/index.html) uses an independently reconstructed exact three-copy partial-transpose certificate: 62 equations, 72 invariant coordinates and 171 positive rational LDL pivots. The [sixth-order proof](research/n17-t6/index.html) is a short rank-one-complement argument, yielding the general theorem `K((r+1)^2+1,2r)=0` for every integer `r>=1`. The existing radical AC4 state therefore has globally maximum order four. Both arguments cover arbitrary complex states.

## Eighteen qubits: fifth and sixth order settled

**K(18,5)=1 and K(18,t)=0 for every t>=6.** The [existing four-component Julia-source witness](research/n18-t5/index.html) is verified exactly, and twelve rectangular tensor images per basis vector give `12k<=18`. The [sixth-order exclusion](research/n18-t6/index.html) is an exact positive identity involving the three-copy trace and total-spin probabilities at J=1,3,4,6. Its independent standard-library verifier reconstructs all 80 equations in 84 coordinates, their four-dimensional kernel, and all physical probabilities. The maximum pure AC order is five. The fourth-order subspace maximum is also one, by the separate auxiliary-spin proof below.

## Nineteen qubits: fifth and sixth order settled

**K(19,t)=0 for every t>=5.** An [exact four-copy matrix-positivity certificate](research/n19-t5/index.html) excludes fifth order for arbitrary complex pure states. Its independent verifier reconstructs 1380 contraction equations and checks 277 positive rational pivots. A separate [three-copy probability identity](research/n19-t6/index.html) excludes sixth order. A new [five-component radical AC4 state](research/n19-t4/index.html) attains the globally maximum pure-state order four; the maximum fourth-order subspace dimension is also one, by the separate exact four-copy plane exclusion below. Both catalogues include the same witness and exclusions.

## Eighteen qubits: fourth-order maximum settled

**K(18,4)=1.** A [short auxiliary-spin-three proof](research/n18-t4/index.html) gives `7k<=13`, excluding every fourth-order AC plane. The existing exact AC5 line attains dimension one. More generally, **7K(N,t)<=N-5 for N>=16 and t>=4**; the builder applies this bound to every applicable entry. The physical spin matrices, exact partial trace, spectrum and attaining state are independently verified. At N=19 this bound alone gives at most two; the independent four-copy certificate below sharpens the maximum to one.


## Nineteen and twenty qubits: all displayed cases settled

The complete rows for orders one through seven are **7, 4, 2, 1, 0, 0, 0** at N=19 and **7, 4, 2, 1, 1, 0, 0** at N=20. The five formerly open cells are now rigorous equalities.

[Exact radical planes at N=19](research/n19-t3/index.html) and [N=20](research/n20-t3/index.html) satisfy all 64 ordered three-qubit reduced-matrix identities. Their whole complex spans are AC3. The new [general auxiliary-spin-five-halves theorem](research/t3-spin5half-inertia/index.html) proves **6K(N,t) ≤ N−4 for every N≥19 and t≥3**. This excludes dimension three at both sizes; the builder applies the theorem to every applicable catalogue entry.

The fourth-order maxima are both one. A hypothetical plane produces a positive normalized moment on the four-copy Young space of shape (2,2). Independently rebuilt physical contraction equations admit a strictly positive rational matrix combination, forcing the moment to vanish. The [N=19 certificate](research/n19-t4/exclusion/index.html) uses 3535 selected equations and 448 positive rational pivots; the [N=20 certificate](research/n20-t4/exclusion/index.html) uses 4253 equations and 518 pivots. Existing exact pure states attain the dimension-one lower bounds. These exclusions cover arbitrary complex subspaces, with no restriction to symmetric constellations or sparse supports.

The [N=20 sixth-order exclusion](research/n20-t6/index.html) uses an independently reconstructed three-copy matrix-positivity identity, with 95 selected equations and 61 positive rational pivots. Thus K(20,t)=0 for all t≥6, and the existing AC5 state has globally maximum pure-state order five. The AC-state catalogue is updated consistently. All proofs, exact coefficients, standalone verifiers, and certified Majorana geometry are bundled.


## Numerical survey at N = 21–24

The [numerical survey](research/n21-24-numerical/index.html) records 3,944 search starts and eleven accepted bases. The largest dimensions found, ordered by t = 1,…,7, are `(8,4,2,1,–,–,–)`, `(8,4,2,1,1,–,–)`, `(8,4,2,1,1,–,–)` and `(9,5,3,1,1,1,1)` for N = 21,22,23,24. A dash means that no candidate was found; it is not a proved zero. First-order dimensions and the previously certified N24 cases retain their exact status.

New numerical examples include dimension four at N21–23, order two, and dimension three at N24, order three. The local Liège source supplies a seventh-order N24 pure-state candidate. Independent 80-digit normalized-Dicke splitting checks all 1,670 ordered scalar reduced entries through the claimed orders, plus every Gram matrix, using only the saved decimal coefficients. Maximum errors are below 1e-14. All components and provenance are bundled; their Majorana representations are available in the catalogue.

The numerical candidates appear as `Num. k` and never raise rigorous lower bounds. Failed higher-dimension and higher-order targets remain open and link to their restart records. For N21 and N23, the existing general first-order construction is instantiated and verified exactly, adding complete rows to the catalogue. The first-order maxima are both eight.


## Majorana barycenters inside maximal first-order AC subspaces

[The N=2–20 survey](research/t1-barycenter-survey/index.html) maximizes the unweighted unit-star barycenter distance within one specified analytical maximal AC1 subspace per size. N=2,3,4 have proved zero maxima; all larger displayed maxima remain numerical candidates. The report includes full-precision state components, interactive Majorana spheres and comparison with the 2014 paper's P_N family. It records 7,039 new completed search starts and reuses 278 previous N12 starts. No AC-dimension grid values change.

Fifteen candidates are defined by isolated simple Wronskian roots, independently certified with exact rational interval arithmetic. The N6 candidate uses an explicit radical-parameter double-root construction. Independent 110-digit checks reconstruct all 19 states and constellations. Global displacement optimality is not proved for N>=5.

Reproduce the exported data and certificates with `research/t1-barycenter-survey/assemble.py`, `verify_independently.py` and `certify_roots.py`; the grouped search scripts and raw logs are included.


## Whole-subspace visualization

The selected representative now has a **View whole subspace** entry. It constructs the orthogonal projector onto the saved Dicke span and uses the trace-one operator `rho_V=P_V/d`, where `d` is the span dimension. Numerical orthonormalization changes only the coordinates used to describe that span; it does not modify the source data or upgrade a numerical AC construction to a proof.

The Husimi view displays `<n^N|rho_V|n^N>` on the physical sphere, with a colored unit-sphere view and a radial option. Its angular mean is `1/(N+1)`, and its operator purity is `1/d`. Camera rotations and axis alignment affect only the display frame.

The operator-constellation view implements the orthonormal spin-tensor decomposition and rank-by-rank Majorana polynomials of E. Serrano-Ensástiga and D. Braun, *Majorana representation for mixed states*, Physical Review A **101**, 022332 (2020), https://doi.org/10.1103/PhysRevA.101.022332 (arXiv:1909.07740, Section III). A nonzero rank `l` has `2l` antipodal stars and physical radius `w_l=||rho_l||`. The signed multipole coefficients and an oriented representative of the subconstellation class preserve the information lost by displaying only unsigned antipodal axes. The scalar rank zero has no stars. A vanished rank has no constellation.

Symmetry axes are tested against the complete operator, including its signed blocks. A symmetry of a single unsigned rank constellation need not be a symmetry of the subspace. The displayed numerical axes and root diagnostics are separate from the catalogue's exact whole-subspace AC proofs. Small low-order multipoles of the projector alone are not a whole-subspace AC certificate; all compressed tensors `P_V T_lm P_V` must vanish at ranks `1,...,t`.
