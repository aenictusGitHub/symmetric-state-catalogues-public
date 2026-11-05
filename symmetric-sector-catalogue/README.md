# Extremal symmetric states

A static research catalogue of pure symmetric multiqubit states extremizing Pauli sector lengths, inspired by [Felix Huber's AME table](https://huberfe.github.io/ame).

The catalogue contains **154 records for 2–12 qubits**, drawn from the 12 supplied PDFs. It includes minimum/maximum tables, evidence filters, a proof-and-gap explorer, exact formulas, representative Dicke states, downloadable records, general constructions and a bibliography with PDF page links.

Every record has a **View Majorana sphere** button for its selected representative. Drag or use the arrow keys to rotate, adjust the zoom, inspect multiplicities, and show symmetry axes and geometric point-group labels. The sphere uses the same convention and controls as the AC and degeneracy catalogues. Downloading a record includes its constellation; all coordinates are also in `data/majorana.json`. The representative panel shows its full point group and proper rotational subgroup without opening the sphere.

The homogeneous polynomial is `P(X,Y) = Σ_r (−1)^r √binom(n,r) c_r X^(n−r)Y^r`; a root `z = Y/X` represents the spinor ratio `b/a = 1/z`, with zero at the south pole and infinity at the north pole. Exact factorization determines repeated points when an exact formula is available. Otherwise the sphere displays roots of the saved decimal coefficients without clustering nearby finite roots. Coordinates and symmetry groups are numerical geometry; the viewer preserves each record's existing optimality status.

## Open locally

Open `index.html` in a web browser. All fonts needed for the equations and all website assets are included locally. The catalogue also works without an internet connection; external DOI links require one.

For a local web address, run this from the website directory:

```sh
python3 -m http.server 8765 --bind 127.0.0.1
```

Then visit `http://127.0.0.1:8765/`. Individual results have links such as `#n=9&k=8&mode=max`.

## Publish on GitHub Pages

1. Upload the **contents of this website directory** to the root of the intended GitHub repository, including `.nojekyll`, `data/`, `vendor/` and `references/`.
2. In the repository, open **Settings → Pages**.
3. Select **Deploy from a branch**, the intended branch (usually `main`), and **/(root)**. Save.

There is no build step, server, API key, package installation or absolute URL configuration. Relative links work under a GitHub project URL. These steps follow the [GitHub Pages publishing documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

`references/` contains copies of all supplied PDFs, including the working manuscript, as the website's downloadable source record. Publication has not been performed.

## Scientific conventions

- `S_k` is the **sum** of squared expectations of all weight-`k` Pauli strings, without binomial normalization. `S_0=1`; a pure state has total sector sum `2^n`.
- States belong to the fully symmetric spin-`n/2` subspace. Dicke amplitudes are ordered by excitation number `r=0,...,n`.
- `proven` means a cited theorem/bound plus an attaining symmetric state, or the short analytic deduction explicitly stated in the entry.
- `reported` means an exact candidate or threshold whose global proof is not reviewable in the supplied material. `numerical` means a numerical candidate. `bound` does not assert attainment.
- Empty cells mean **not catalogued**, not nonexistence or a verified open problem.
- A representative is not a classification of all extremizing states. Collective rotations preserve every sector length.

The working draft distinguishes conjectural entries in its own Appendix A. The research update located the pinned public threshold supplement and added independent exact certificates. It proves 43 further records: all 154 are now analytically established. Details and reproducible checks are in [research/README.md](research/README.md) and the [proof page](research/index.html). Remaining limitations are retained in `data/catalogue.json` under `sourceCaveats`. The nine-qubit eighth-sector theorem is attributed to the working draft; its separate sum-of-squares certificate repository was not rerun. The local October 2026 anticoherence manuscript is distinguished from its related public arXiv version.

## Data and maintenance

`data/catalogue.json` is the canonical editable dataset. Each record includes `(n,k,mode)`, exact and numeric values, evidence status, state formula, available coefficients, notes and source references. Source metadata includes SHA-256 hashes of the supplied PDFs. `data/catalogue.csv` is a flat summary; the JSON retains full coefficients and references.

After changing the JSON, regenerate the browser data and CSV:

```sh
python3 scripts/build_data.py
```

The generated `data/catalogue.js` lets the site open directly from a local file without a web server. No browser data is stored or transmitted.

## Verification

For the optional scientific check, use Python 3.10+ with NumPy:

```sh
python3 scripts/verify_states.py --self-test data/catalogue.json
```

The checker forms reduced density matrices from normalized Dicke amplitudes and obtains sector lengths by binomial inversion. Self-tests compare against direct Pauli-string enumeration for complex symmetric states and analytic coherent/GHZ identities.

On the supplied dataset, **154 representative checks pass, none are skipped and none fail**. Entries use their declared numerical tolerance; exact algebraic reconstructions use `2e-7` for their floating-point coefficient evaluations. These are checks of the representative values, not proofs of global optimality. Every record now has an explicit representative. The research update additionally verifies 154 rational bound identities and 27 exact witness cases, plus three fourth-sector constructions.

The website was also checked for mathematical rendering, source links, minimum/maximum selection, evidence filtering, missing-result handling and responsive layout.

After changing representative data, regenerate the Majorana coordinates with `python3 scripts/build_majorana.py` (SymPy and NumPy). Independent reconstruction and viewer checks are included:

```sh
python3 research/majorana-viewer/verify_coordinates.py
node research/majorana-viewer/verify_viewer.cjs
node research/majorana-viewer/verify_app.cjs
```

The coordinate verifier uses mpmath; the JavaScript checks use the locally bundled KaTeX. Saved reports are in `research/majorana-viewer/`.

## Included renderer

KaTeX 0.16.22 is included in `vendor/katex/`; its MIT license is preserved there. The website does not require a CDN. Scientific sources retain their respective attribution.

The follow-up reconstructs exact algebraic stationary states for the S6 minimum candidates at n=10,11,12 and an explicit radical state for the n=9 S5 maximum candidate. See [algebraic S6 states](research/algebraic-s6/index.html) and [the radical S5 state](research/algebraic-s5-n9/index.html). Global optimality has since been proved for all four, including the eleven-qubit S6 minimum. The three defining parameter polynomials have degrees14,11,12; the n=9 state uses nested radicals over Q(sqrt(61)).


The ten-qubit S6 entry now exposes its **exact attained value** next to the decimal.
Its value polynomial is irreducible of degree 14, with Galois group S14; therefore
there is no expression by radicals over the rationals. The exact value and state,
root-isolation certificate, and finite-field proof are in the
[analytical appendix](research/algebraic-s6/index.html#n10-value).
Reproduce the exact arithmetic with
`python3 research/algebraic-s6/n10/verify_value.py` (SymPy).
This strengthens the description of the candidate, not its global-optimality status.


The [six-qubit fourth-sector minimum](research/n6-s4-minimum/index.html) is now
proved to equal **45/7** by an exact rational degree-eight sum-of-squares certificate.
The proof covers arbitrary complex symmetric states and includes two inequivalent
equality states. The global and equality-condition verifiers require only ordinary
Python 3; their full coefficient and positivity checks use exact rational arithmetic.


The [seven-qubit sixth-sector minimum](research/n7-s6-minimum/index.html) is
proved to equal **7** by an exact rational degree-six sum-of-squares certificate.
Every minimizer is a collective rotation of a normalized generalized GHZ state.
The proof also gives the sharp tradeoff S2+S4 <= 56 and an equivalent inequality
between reduced-state purities. The global, even-sector and complete
equality-classification checks use only the Python standard library.


## Seven-qubit fourth-sector maximum (3 November 2026)

The global maximum is exactly **265/7**. An exact rational degree-six
sum-of-squares certificate proves the bound for arbitrary complex symmetric
states. Six positive definite Gram matrices yield 60 weighted absolute squares.
The existing sparse Dicke representative attains the bound; so does its whole
complex span with its time reversal, after normalization. Completeness of this
attaining family is not claimed.

[Read the proof and download the exact certificate](research/n7-s4-maximum/index.html).

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

[Read the proof and download the certificate](research/n7-s5-maximum/index.html).

```sh
python3 research/n7-s5-maximum/verify_certificate.py
```

The standard-library verifier checks the complete complex-polynomial identity,
exact matrix positivity, W-state attainment and its sector vector, and the
reduced-purity identity, using rational arithmetic throughout.


## All eight-qubit cases settled (3 November 2026)

Exact rational sum-of-squares certificates prove the remaining three extrema:

- [max S4 = 70](research/n8-s4-maximum/index.html): degree 6, 14 positive Gram matrices, 182 weighted absolute squares.
- [max S5 = 147/2](research/n8-s5-maximum/index.html): degree 6, 19 positive Gram matrices, 345 weighted absolute squares.
- [min S6 = 177/14](research/n8-s6-minimum/index.html): degree 8, 41 positive Gram matrices, 3189 weighted absolute squares.

The representatives are D8^0, D8^1, and D8^3, respectively. Each proof folder
contains a standard-library `verify_certificate.py` that independently checks
the full complex-polynomial identity, exact positivity, and Dicke-state attainment.
All eight-qubit entries are now proved. No full equality classification is claimed.

Convexity and symmetric reduction extend max S4 = 70 at n=8 to
max S4 = binomial(n,4) for every n >= 8. The proved n=7 maximum is larger
than binomial(7,4), so the exact eventual threshold is eight.


## Nine-qubit exact rational proofs (3 November 2026)

Four further global extrema are proved for arbitrary complex symmetric states:

- [min S4 = 2754/91](research/n9-s4-minimum/index.html): degree 8, 24 positive Gram matrices, 1575 weighted absolute squares.
- [max S6 = 1204/9](research/n9-s6-maximum/index.html): degree 6, 15 positive Gram matrices, 256 weighted absolute squares.
- [min S6 = 420/11](research/n9-s6-minimum/index.html): degree 8, 30 positive Gram matrices, 1894 weighted absolute squares.
- [max S7 = 1234/9](research/n9-s7-maximum/index.html): degree 6, 19 positive Gram matrices, 438 weighted absolute squares.

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

[Full algebraic proof and reproducible certificate](research/n9-s5-maximum/index.html).

All five previously open nine-qubit cases are now proved. Every catalogue
entry through nine qubits is settled; equality classifications are not implied.


## Ten-qubit exact rational proofs (3 November 2026)

Three further global extrema are proved for arbitrary complex symmetric states:

- [max S6 = 270](research/n10-s6-maximum/index.html): degree 6, 19 positive Gram matrices, 550 weighted absolute squares.
- [min S8 = 21](research/n10-s8-minimum/index.html): degree 8, 42 positive Gram matrices, 5853 weighted absolute squares.
- [max S9 = 880/3](research/n10-s9-maximum/index.html): degree 8, 41 positive Gram matrices, 5080 weighted absolute squares.

Each proof folder contains a standard-library `verify_certificate.py` that
reconstructs the target from reduced density matrices, checks every complex
polynomial coefficient and exact matrix positivity, and verifies the attaining
state. Collective rotations of the representatives also attain the bounds.
No complete equality classification is claimed.


The sixth-sector minimum previously had a small remaining gap. An
[exact degree-eight certificate](research/n10-s6-lower-bound/index.html) improves its global
lower bound from 62 to **58437/800 = 73.04625**. The exactly attained
algebraic candidate is 73.0462518114328286…, leaving a gap below 0.000001812.
This historical bound is superseded by the sharp algebraic proof below.


## Ten-qubit sixth-sector minimum proved (3 November 2026)

The [sharp algebraic proof](research/n10-s6-minimum/index.html) establishes
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

- [max S6 = 462](research/n11-s6-maximum/index.html): degree 6, 13 positive Gram matrices, 200 weighted squares.
- [max S7 = 250060/517](research/n11-s7-maximum/index.html): degree 6, 14 positive Gram matrices, 360 weighted squares.
- [max S8 = 5271/11](research/n11-s8-maximum/index.html): degree 6, 20 positive Gram matrices, 472 weighted squares.
- [min S8 = 77](research/n11-s8-minimum/index.html): degree 6, 13 positive Gram matrices, 208 weighted squares.
- [max S9 = 134645/239](research/n11-s9-maximum/index.html): degree 8, 45 positive Gram matrices, 7891 weighted squares.
- [max S10 = 19459/35](research/n11-s10-maximum/index.html): degree 8, 39 positive Gram matrices, 4622 weighted squares.
- [min S10 = 11](research/n11-s10-minimum/index.html): degree 6, 22 positive Gram matrices, 570 weighted squares.

The [sixth-sector threshold proof](research/n11-s6-maximum/index.html#threshold) establishes **max S6 = binom(N,6) for every N >= 11**. The least eventual threshold is exactly eleven, since max S6 = 270 > 210 at N=10. This also proves the twelve-qubit maximum 924.

The last eleven-qubit gap is now closed: [min S6 = 113.109043566346690991…](research/n11-s6-minimum/index.html) is the exact degree-eleven algebraic global minimum. A degree-ten SOS certificate averaged over Kramers mixing is independently verified using 286 evaluation equations and four exact vanishing constraints in the complete 290-dimensional invariant space. All 22 eleven-qubit extrema, and all 154 catalogue extrema, are proved. The exact minimizing family includes every normalized combination of a witness and its time-reversed partner. Earlier rational lower bounds remain available as historical certificates.


## Twelve-qubit results (4 November 2026)

- [max S8 = 15405/11](research/n12-s8-maximum/index.html): exact global SOS proof and analytical attainment.
- [min S10 = 594/19](research/n12-s10-minimum/index.html): exact global SOS proof and analytical attainment.
- [max S11 = 1270976/1191](research/n12-s11-maximum/index.html): exact global SOS proof and analytical attainment.

All twelve-qubit extrema are now proved. The two algebraic minima are:

- [min S6 = 142.687882135894087969…](research/n12-s6-minimum/index.html), algebraic of degree twelve.
- [min S8 = 223.307697897758835463…](research/n12-s8-minimum/index.html), algebraic of degree fourteen.

Both proofs use exact degree-eight sum-of-squares identities, with independent verification of the full invariant polynomial identity and rigorous positivity. The older rational lower bounds are retained as historical certificates.

The [S8 minimum representative](research/algebraic-s8-n12/index.html) is now analytical: its coefficients and attained value are defined by isolated degree-fourteen algebraic roots. Exact stationarity and global optimality are verified. The main catalogue omits certificate-size statistics; the linked proof pages retain those details.

The eleven-qubit verification reports now display all eleven sector lengths; the previously proved extremal values are unchanged.


## Twelve-qubit ninth-sector maximum settled (4 November 2026)

The exact global maximum is **1280**, attained by `(D12(0) + sqrt(2) D12(6) - D12(12))/2`. The [degree-eight rational SOS proof](research/n12-s9-maximum/index.html) applies to arbitrary complex pure symmetric states. An additional explicit one-parameter family attains the same maximum. The standalone verifier checks the full polynomial identity, positive rational matrices, and exact attainment.
