# Research updates: 28 October – 4 November 2026

This update adds 13 missing maximum records, replaces seven lower-bound-only records with explicit attained candidate values, supplies the four missing witnesses, and proves 43 further records, including the subsequent six- and seven-qubit SOS proofs. All 154 catalogue records are now analytically established. The final proof is the averaged degree-ten certificate for the eleven-qubit sixth-sector minimum.

## Initial global proofs (28 October)

| n | Sector | Direction | Extremal value |
|---|---|---|---|
| 8 | 4 | min | 310/21 |
| 8 | 6 | max | 364/3 |
| 8 | 7 | max | 592/7 |
| 10 | 4 | min | 42 |
| 10 | 5 | max | 252 |
| 10 | 7 | max | 320 |
| 10 | 8 | max | 405 |
| 11 | 4 | min | 66 |
| 11 | 5 | max | 462 |
| 12 | 4 | min | 99 |
| 12 | 7 | max | 946 |
| 12 | 10 | max | 1606 |

The rational certificates also give additional proofs of the already-established n=7 max S6=49 and n=12 max S5=792. `PROOFS.md` derives the two-copy spectral formulas and prints the sharp identities. The only inputs to the sharp two-copy certificates are nonnegative total-spin probabilities and nonnegative squared multipole lengths. No numerical assumption enters these proofs.

The fourth-sector constructions are traced to the [public threshold supplement at commit a992fde3e73455f81445c872a442c640ddb6a11b](https://github.com/aenictusGitHub/n9-s8-exact-certificates/tree/a992fde3e73455f81445c872a442c640ddb6a11b/thresholds). Their normalization and attainment were checked independently. Source checksum and verification metadata are in `threshold-verification.json`; original third-party scripts are not included.

## Verification

Run from the website directory:

```sh
python3 research/certificates/verify_certificates.py
python3 research/verify_exact_witnesses.py
python3 research/verify_minimum_s4.py
python3 scripts/verify_states.py --self-test data/catalogue.json
```

- The first command needs only the standard library. It checks all 154 rational lower/upper bound identities, the nonnegativity of certificate coefficients, and the feasibility of the reported relaxation vertices. A relaxation vertex is not assumed physical.
- The second needs SymPy and verifies 27 exact representative/sector pairs with exact normalization and exact reduced purities. These include all eleven sharp spectral-certificate cases and all 19 exact witnesses recovered in the search, with three overlaps.
- The third independently verifies the three exact n=10,11,12 fourth-sector constructions using SymPy.
- The last uses NumPy and checks all 154 current representatives. Original rounded source states keep their declared rounding tolerance; newly found witnesses have full-precision coefficients.

An exact witness check proves attainment only. A global inequality plus a matching witness proves optimality. Numerical discovery, exact attainment, and a global proof are separate evidence categories throughout the website.

## Reproduce the numerical searches

The search files require NumPy and SciPy, with SymPy for `exact.py`:

```sh
python3 research/search/optimize.py check
python3 research/search/optimize.py 11 6 min 1024
python3 research/search/sparse.py 11 6 min 4
python3 research/search/verify_independent.py
```

The multistart search uses normalized complex Dicke amplitudes, an analytic gradient, and L-BFGS-B. It seeds from the original catalogue (`catalogue-before-update.json`) plus perturbed Dicke/GHZ states and seeded random states. Run files retain the seed, number of starts, termination status, gradient norm and candidate coefficients. Some searches terminate at floating-point line-search precision; independent value checks and subsequent stationary refinement are therefore recorded separately. Difficult cases were revisited with more than 1,000 starts. This is evidence for a candidate, not an exhaustive search or proof.

`verify_independent.py` evaluates Pauli expectations in the full computational Hilbert space, grouping strings only by the symmetric counts of X, Y and Z. This is independent of the reduced-purity objective used in optimization. `refine_shapes.py` refines the simple real-support patterns of the three difficult minima, using saved support-search results.

`certificates/lp_bounds.py` reproduces rational-certificate discovery with SciPy and SymPy. The discovery optimizer is not part of the exact proof verifier.

## Remaining open proof gaps in this catalogue

The 3 candidate optima have separate certified lower or upper bounds. They are listed on `index.html` and in the canonical catalogue JSON. The remaining cases are:

- n=11 minimum S6: candidate 113.109043566347; certified lower bound 565539/5000.
- n=12 minimum S6: candidate 142.687882135894; certified lower bound 713439361/5000000.
- n=12 minimum S8: candidate 223.307697897759; certified lower bound 22330768799/100000000.

“Remaining” refers to the proofs established in this catalogue, not a claim that no proof exists elsewhere. Bounded three-copy spectral and multiplier searches did not produce additional exact certificates, and no result was upgraded on that basis.

All convention choices and source references remain in `data/catalogue.json`. The initial dataset is preserved as `catalogue-before-update.json` for comparison and reproducibility.


## Six-qubit fourth-sector minimum proved (3 November 2026)

The global minimum is exactly **45/7**. The degree-eight rational sum-of-squares
identity applies to arbitrary complex symmetric states. Independent exact checks
verify all polynomial coefficients and all 18 positive Gram matrices. Two
rotation-inequivalent equality states are supplied, along with necessary and
sufficient spin-moment conditions for equality. See the
[proof and explicit certificate](n6-s4-minimum/index.html).

```sh
python3 research/n6-s4-minimum/verify_certificate.py
python3 research/n6-s4-minimum/verify_equality_conditions.py
```

Both verifiers use only the Python standard library. The old two-copy certificate
27/5 remains valid in `certificates/lp-bounds.json`; the new nonlinear certificate
is stored separately and brought the catalogue to 124 proved entries and 30
remaining candidates before the seven-qubit result below.


## Seven-qubit sixth-sector minimum and all equality states (3 November 2026)

The global minimum is exactly **7**, with equality precisely for collective
rotations of normalized generalized GHZ states. An exact rational degree-six
sum-of-squares certificate and an independent complete equality classification
settle this case. Sharp corollaries are S2+S4 <= 56 and the reduced-purity
inequality 13 Tr(rho1²)+20 Tr(rho3²) >= 33 Tr(rho2²).

[Read the full proof and explicit certificate](n7-s6-minimum/index.html).

```sh
python3 research/n7-s6-minimum/verify_certificate.py
python3 research/n7-s6-minimum/verify_even_sectors.py
python3 research/n7-s6-minimum/verify_classification.py
```

These three checks use only the Python standard library and exact fractions.
The old 20/3 linear-relaxation certificate remains valid and is kept separately.
The catalogue now contains 151 proved records and 3 remaining candidates.


## Seven-qubit fourth-sector maximum (3 November 2026)

The global maximum is exactly **265/7**. An exact rational degree-six
sum-of-squares certificate proves the bound for arbitrary complex symmetric
states. Six positive definite Gram matrices yield 60 weighted absolute squares.
The existing sparse Dicke representative attains the bound; so does its whole
complex span with its time reversal, after normalization. Completeness of this
attaining family is not claimed.

[Read the proof and download the exact certificate](n7-s4-maximum/index.html).

```sh
python3 research/n7-s4-maximum/verify_certificate.py
```

The verifier uses only the Python standard library and exact fractions.
It checks the complete polynomial identity, positivity, the exact attaining
state and its sector vector, and the stated continuous attaining family.


## Seven-qubit fifth-sector maximum (3 November 2026)

The global maximum is exactly **267/7**, attained by the W state D7^1 and its
collective rotations. A rational degree-six sum-of-squares certificate proves
the bound for arbitrary complex symmetric states. Fifteen positive definite
Gram matrices of size at most three give 226 weighted absolute squares.

An equivalent sharp purity inequality is
5 Tr(rho1²) - 4 Tr(rho2²) <= 69/49. The exact W-state purities are 37/49 and 29/49.
Completeness of the attaining orbit is not claimed.

[Read the proof and download the certificate](n7-s5-maximum/index.html).

```sh
python3 research/n7-s5-maximum/verify_certificate.py
```

The standard-library verifier checks the complete complex-polynomial identity,
exact matrix positivity, W-state attainment and its sector vector, and the
reduced-purity identity, using rational arithmetic throughout.


## All eight-qubit cases settled (3 November 2026)

Exact rational sum-of-squares certificates prove the remaining three extrema:

- [max S4 = 70](n8-s4-maximum/index.html): degree 6, 14 positive Gram matrices, 182 weighted absolute squares.
- [max S5 = 147/2](n8-s5-maximum/index.html): degree 6, 19 positive Gram matrices, 345 weighted absolute squares.
- [min S6 = 177/14](n8-s6-minimum/index.html): degree 8, 41 positive Gram matrices, 3189 weighted absolute squares.

The representatives are D8^0, D8^1, and D8^3, respectively. Each proof folder
contains a standard-library `verify_certificate.py` that independently checks
the full complex-polynomial identity, exact positivity, and Dicke-state attainment.
All eight-qubit entries are now proved. No full equality classification is claimed.

Convexity and symmetric reduction extend max S4 = 70 at n=8 to
max S4 = binomial(n,4) for every n >= 8. The proved n=7 maximum is larger
than binomial(7,4), so the exact eventual threshold is eight.


## Nine-qubit exact rational proofs (3 November 2026)

Four further global extrema are proved for arbitrary complex symmetric states:

- [min S4 = 2754/91](n9-s4-minimum/index.html): degree 8, 24 positive Gram matrices, 1575 weighted absolute squares.
- [max S6 = 1204/9](n9-s6-maximum/index.html): degree 6, 15 positive Gram matrices, 256 weighted absolute squares.
- [min S6 = 420/11](n9-s6-minimum/index.html): degree 8, 30 positive Gram matrices, 1894 weighted absolute squares.
- [max S7 = 1234/9](n9-s7-maximum/index.html): degree 6, 19 positive Gram matrices, 438 weighted absolute squares.

Each proof folder contains a standard-library `verify_certificate.py` that
reconstructs the target from reduced density matrices, checks every complex
polynomial coefficient and exact matrix positivity, and verifies the attaining
state. Even-sector certificates additionally verify the full Kramers span.
No complete equality classification is claimed.


## Nine-qubit fifth-sector maximum: all n <= 9 entries settled

The global S5 maximum is exactly **(130928 + 6832 sqrt(61))/1461**,
attained by the explicit nested-radical state on Dicke support 0, 3, 6, 9.
An exact degree-eight SOS certificate over Q(sqrt(61)) proves the bound for
arbitrary complex symmetric states. The independent standard-library verifier
checks every coefficient, every Gram matrix, and exact radical attainment.

[Full algebraic proof and reproducible certificate](n9-s5-maximum/index.html).

All five previously open nine-qubit cases are now proved. Every catalogue
entry through nine qubits is settled; equality classifications are not implied.


## Ten-qubit exact rational proofs (3 November 2026)

Three further global extrema are proved for arbitrary complex symmetric states:

- [max S6 = 270](n10-s6-maximum/index.html): degree 6, 19 positive Gram matrices, 550 weighted absolute squares.
- [min S8 = 21](n10-s8-minimum/index.html): degree 8, 42 positive Gram matrices, 5853 weighted absolute squares.
- [max S9 = 880/3](n10-s9-maximum/index.html): degree 8, 41 positive Gram matrices, 5080 weighted absolute squares.

Each proof folder contains a standard-library `verify_certificate.py` that
reconstructs the target from reduced density matrices, checks every complex
polynomial coefficient and exact matrix positivity, and verifies the attaining
state. Collective rotations of the representatives also attain the bounds.
No complete equality classification is claimed.


The sixth-sector minimum previously had a small remaining gap. An
[exact degree-eight certificate](n10-s6-lower-bound/index.html) improves its global
lower bound from 62 to **58437/800 = 73.04625**. The exactly attained
algebraic candidate is 73.0462518114328286…, leaving a gap below 0.000001812.
This historical bound is superseded by the sharp algebraic proof below.


## Ten-qubit sixth-sector minimum proved (3 November 2026)

The [sharp algebraic proof](n10-s6-minimum/index.html) establishes
min S6 = s_star = 73.0462518114328286265592652767668… exactly, where s_star
is the previously reconstructed degree-14 algebraic number. All twenty
ten-qubit sector extrema are now proved.

The degree-eight SOS certificate has 34 positive Gram blocks. Exact covariance
and highest-weight checks establish a complete 130-dimensional invariant space;
126 coefficient equations and four exact vanishing constraints prove the full
identity. Directed-rounding Arb arithmetic certifies injectivity, the algebraic
linear solve and all positive LDL pivots, with no numerical tolerance assumption.
The verifier requires python-flint; no optimizer, NumPy or SciPy is needed.


## Eleven-qubit exact bounds (4 November 2026)

Seven previously open global extrema are proved for arbitrary complex symmetric states:

- [max S6 = 462](n11-s6-maximum/index.html): degree 6, 13 positive Gram matrices, 200 weighted squares.
- [max S7 = 250060/517](n11-s7-maximum/index.html): degree 6, 14 positive Gram matrices, 360 weighted squares.
- [max S8 = 5271/11](n11-s8-maximum/index.html): degree 6, 20 positive Gram matrices, 472 weighted squares.
- [min S8 = 77](n11-s8-minimum/index.html): degree 6, 13 positive Gram matrices, 208 weighted squares.
- [max S9 = 134645/239](n11-s9-maximum/index.html): degree 8, 45 positive Gram matrices, 7891 weighted squares.
- [max S10 = 19459/35](n11-s10-maximum/index.html): degree 8, 39 positive Gram matrices, 4622 weighted squares.
- [min S10 = 11](n11-s10-minimum/index.html): degree 6, 22 positive Gram matrices, 570 weighted squares.

The [sixth-sector threshold proof](n11-s6-maximum/index.html#threshold) establishes **max S6 = binom(N,6) for every N >= 11**. The least eventual threshold is exactly eleven, since max S6 = 270 > 210 at N=10. This also proves the twelve-qubit maximum 924.

The last eleven-qubit gap is now closed: [min S6 = 113.109043566346690991…](n11-s6-minimum/index.html) is the exact degree-eleven algebraic global minimum. A degree-ten SOS certificate averaged over Kramers mixing is independently verified using 286 evaluation equations and four exact vanishing constraints in the complete 290-dimensional invariant space. All 22 eleven-qubit extrema, and all 154 catalogue extrema, are proved. The exact minimizing family includes every normalized combination of a witness and its time-reversed partner. Earlier rational lower bounds remain available as historical certificates.


## Twelve-qubit results (4 November 2026)

- [max S8 = 15405/11](n12-s8-maximum/index.html): exact global SOS proof and analytical attainment.
- [min S10 = 594/19](n12-s10-minimum/index.html): exact global SOS proof and analytical attainment.
- [max S11 = 1270976/1191](n12-s11-maximum/index.html): exact global SOS proof and analytical attainment.

All twelve-qubit extrema are now proved. The two algebraic minima are:

- [min S6 = 142.687882135894087969…](n12-s6-minimum/index.html), algebraic of degree twelve.
- [min S8 = 223.307697897758835463…](n12-s8-minimum/index.html), algebraic of degree fourteen.

Both proofs use exact degree-eight sum-of-squares identities, with independent verification of the full invariant polynomial identity and rigorous positivity. The older rational lower bounds are retained as historical certificates.

The [S8 minimum representative](algebraic-s8-n12/index.html) is now analytical: its coefficients and attained value are defined by isolated degree-fourteen algebraic roots. Exact stationarity and global optimality are verified. The main catalogue omits certificate-size statistics; the linked proof pages retain those details.

The eleven-qubit verification reports now display all eleven sector lengths; the previously proved extremal values are unchanged.


## Twelve-qubit ninth-sector maximum settled (4 November 2026)

The exact global maximum is **1280**, attained by `(D12(0) + sqrt(2) D12(6) - D12(12))/2`. The [degree-eight rational SOS proof](n12-s9-maximum/index.html) applies to arbitrary complex pure symmetric states. An additional explicit one-parameter family attains the same maximum. The standalone verifier checks the full polynomial identity, positive rational matrices, and exact attainment.
