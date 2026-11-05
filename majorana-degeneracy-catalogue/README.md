# Anticoherence by full Majorana multiplicity partition

Open `index.html` in a browser. This static catalogue uses bundled formulas, references and data, and needs no server or build step for normal use.

Each entry is indexed by the **qubit number n and the complete Majorana multiplicity partition λ**. Its value is the **highest achievable anticoherence order**, maximized over all pure symmetric states with exactly that partition. Select a qubit number to browse its partitions, or choose all sizes.

Choose **View Majorana sphere** above the table or beside a representative to open the interactive constellation. Drag to rotate, use the arrow buttons or keyboard, adjust zoom, and switch representatives inside the viewer. All 197 attached states have coordinates; partitions without an explicit representative display an explanation instead of an empty sphere.

The **Symmetry axes** overlay shows proper rotational symmetries of the weighted star set, preserving every point multiplicity. Each unoriented axis carries its highest order, such as C₂ or C₅; the order selector filters the drawing. For the ten-qubit cyclic representative, there is one fivefold axis and five twofold axes. For dense polyhedral axis sets, the viewer initially shows the highest-order axes; choose “All orders” to see the rest. Coordinate axes have a separate toggle. An antipodal balanced constellation shows its continuous rotation axis and the perpendicular plane of twofold axes.

The table, representative card and sphere viewer also show the **full point group** and its **proper rotational subgroup**. For example, the icosahedral representative has point group Iₕ with 120 symmetries and rotational subgroup I with 60 rotations. Both groups preserve point multiplicities: the weighted double tetrahedron with its simple dual has T_d / T, even though its unweighted positions form a cube. The table lists the groups of all attached representatives; a partition need not have a unique group. Rear stars have a warm fill, darker labels and dashed outlines.

Symmetry detection is numerical, using a unit-sphere distance tolerance of 10⁻⁶. It does not change the exact AC claims or certify algebraic symmetry. The offline modules are `majorana-symmetry.js`, `symmetry-groups.js` and `point-groups.js`. Geometric, closure, improper-coset and fixture checks are recorded in `research/symmetry-validation.json`, `research/symmetry-group-validation.json` and `research/point-group-validation.json`. [Download the per-state groups](data/symmetry-groups.json); regenerate them with `node scripts/build_symmetry_groups.cjs` after changing the state data.

The main table follows the requested full-partition organization. The earlier n × AC-order table of maximum single-point multiplicity remains in the downloadable data as `legacyEntries`, supplying some universal upper bounds.

## Coverage

The catalogue contains **387 full partitions**: all **270 integer partitions for 2 ≤ n ≤ 12**, plus **117 source-observed partitions for 13 ≤ n ≤ 100**. The larger sizes are selected observations, not complete partition enumerations.

The complete small-size core contains:

- **126 proved positive maxima**.
- **112 proved zeros**: no state with that partition has any positive AC order.
- **32 rigorous intervals**: the highest achievable order remains undetermined within the displayed lower and upper bounds.

There are **48 exact representatives** and **144 numerical representatives**. The exact states cover 40 full partitions: 41 have radical Dicke coefficients, six use isolated roots of univariate integer polynomials, and one eleven-qubit state is specified by a rational polynomial system and a certified isolating box. Other rigorous positive lower bounds follow from the finite-filtering existence theorem; no coefficient vector is claimed where none is attached.

The numerical inputs correspond to **145 records listed by the Liège libraries**, of which **144 have available coordinate files**. The package preserves **243 original files**: 144 coordinate files and 99 Dicke files. The listed `acs.12.2.3.txt` coordinate file is missing. Files supplied locally were matched to the [general library](https://www.oq.ulg.ac.be/libspinACstates.html) and [degenerate-state library](https://www.oq.ulg.ac.be/libspinACstatesdeg.html) inventories; current online byte identity could not be checked while the host was unreachable. Provenance, hashes and the missing record are recorded in `research/liege-library/manifest.json`.

Numerical observations do not strengthen a rigorous lower bound or establish a maximum. Their reported orders and residuals are shown separately from exact results.

## What the values mean

A partition such as `(2²,1³)` means two double points and three distinct simple points, for a total of seven qubits. It does not specify their locations. The maximum is taken over all constellations with that exact partition, not over every pattern sharing the same largest part.

A single number is a proved maximum. An interval `a–b` means a state exists at order a and no state in the partition can exceed b. It does not assert attainability at b. The value zero uses the harmless order-zero convention: symmetric states with the partition exist, but none is first-order AC. Positive-order nonexistence propagates to all higher orders.

For example:

| n | Full partition | Highest AC order |
|---|---|---|
| 7 | (2,1⁵) | 2, proved |
| 7 | (2²,1³) | 1, proved |
| 12 | (2⁶) | 3, proved |
| 12 | (3⁴) | 2, proved |
| 12 | (2⁴,1⁴) | between 2 and 3 |
| 12 | (1¹²) | 5, proved |

The selected state always matches the entire partition. Some partitions have comparison representatives of lower AC order; their verified order is displayed separately. Failure of one representative at the next order never serves as a nonexistence proof for its partition.

## Improved AC orders

Nine partitions now have stronger rigorous lower bounds. All new witnesses are exactly second-order AC, with their full multiplicities certified.

| n | Full partition | Previous bounds | Updated bounds |
|---|---|---|---|
| 8 | (2,1⁶) | 1–2 | **2, proved** |
| 9 | (2³,1³) | 1–2 | **2, proved** |
| 10 | (2,1⁸) | 1–2 | **2, proved** |
| 10 | (2³,1⁴) | 1–2 | **2, proved** |
| 11 | (2²,1⁷) | 1–3 | **2–3** |
| 11 | (2³,1⁵) | 1–3 | **2–3** |
| 11 | (2⁴,1³) | 1–3 | **2–3** |
| 12 | (3,2³,1³) | 1–2 | **2, proved** |
| 12 | (3³,1³) | 1–2 | **2, proved** |

Three sparse states use radical coefficients. Six cyclic constructions use factored Majorana polynomials and isolated algebraic parameters. Their certificates check the original moment equations exactly, positivity of normalization, and distinctness of the polynomial factors. Rounded displayed coefficients are not the proof. See the [sparse certificates](research/ac-order-sparse.json), [cyclic and two-pole certificates](research/ac-order-constructions.json), and [independent audit](research/ac-order-improvements-audit.json).

For the selected twelve-qubit partition `(2⁴,1⁴)`, the rigorous range remains **2–3**. A focused follow-up now proves exact exclusions for four doubled points in two antipodal pairs, and for proper rotational symmetry groups other than C₁ or C₂. These leave arbitrary asymmetric or solely twofold-symmetric geometries unresolved. The [readable proof note](research/n12-four-double/index.html) gives the arguments and independent exact verifiers.

The follow-up also constructs an exact sequence inside `(2⁴,1⁴)` tending to a third-order state of partition `(2⁶)`. Therefore the infimum of the third-order defect is zero even within the target partition, without proving that it is attained. A further 2,866 complex-coordinate searches/root solves and 1,500 antipodal searches found no validated third-order witness. See the [search summary](research/n12-four-double/search-summary.json), [antipodal search](research/n12-four-double/search_antipodal.json), and [earlier search record](research/high-order-search.json). Numerical failures are not nonexistence proofs.

## Exact higher-size witnesses

Five analytic formulas from the supplied `ACLiege.jl` now strengthen the all-simple partitions: the rigorous range becomes **3–5 at n=13**, **3–6 at n=15**, and **5–7 at n=16,18,20**. Exact reduced matrices certify the AC orders, and algebraic square-free factorization certifies every Majorana point is simple. The existing 144 numerical representatives remain separate and unchanged. These are lower-bound improvements, not new global maximality proofs. See [exact formulas](research/acliege/exact-states.json), [verification](research/acliege/verification.json), and the companion [AC-subspace catalogue](../ac-subspace-catalogue/index.html).

## Refined numerical components

All **144 primary numerical representatives** now satisfy a Hilbert–Schmidt residual below **10⁻⁶⁰ at their original reported orders**, checked after reading back the exported components. The largest residual is approximately **5.872 × 10⁻⁶⁶**. **135 representatives were improved; nine exact balanced Dicke representatives were already accurate and unchanged.** That refinement kept the source orders fixed. The separate exact constructions described below improve rigorous bounds in nine small-size partitions.

Original files remain intact. The downloads contain 85 significant digits per real and imaginary component. Repeated-pole structural zeros are preserved. Five singular second-order cases use their forced two-term family, removing tiny forbidden amplitudes while retaining the original nonzero phases. The largest component change is about 1.136 × 10⁻⁷.

The **16 unmatched original Dicke counterparts** are also refined, as separate download-only alternatives. They have no assigned partition or sphere plot and do not replace the coordinate-consistent primary representatives. The other 83 original Dicke vectors matched their coordinate counterparts to high precision after convention changes and a global phase.

The displayed refined sphere points match the refined primary components. Independent root-count checks certify the full partitions of the **saved numerical vectors**, including 5,437 simple-root neighborhoods and 238 fixed pole roots. These checks do **not** certify exact anticoherence. Printed precision and small residuals remain numerical statements.

[Download all refined components](research/liege-refinement/refined-components.zip) · [Methods and reproducible files](research/liege-refinement/README.md) · [Coefficient download manifest](research/liege-refinement/download-manifest.json)

## Proof methods

The first-order classification is complete. If the largest part is less than n/2, a finite invertible collective filter produces a first-order AC state without changing any multiplicity. The norm-minimization proof establishes existence directly. If the largest part exceeds n/2, first order is impossible. At equality, only the balanced two-point partition `(n/2,n/2)` permits first order, attained by the balanced Dicke state.

Higher-order upper bounds combine the exact AC nonexistence results from the companion catalogue with moment bounds on point multiplicities. At n=3g+1, saturation of the second-order bound forces the full partition `(g,1^(n−g))`, giving a restriction beyond the largest part alone. Exact representatives establish stronger lower bounds in their own partitions.

The proof section identifies independently derived arguments and provides the source context. Numerical source surveys are not used as proofs. The new repeated-tetrahedron witnesses at twelve qubits have independently checked reduced matrices and exact polynomial factorizations. The eleven-qubit double-point state has exact existence and multiplicity certificates with independent audits.

## Conventions

The normalized Dicke coefficient c_r has r excitations. The homogeneous Majorana polynomial is

```math
P(X,Y)=\sum_{r=0}^{n}(-1)^r\sqrt{\binom nr}\,c_rX^{n-r}Y^r.
```

For affine coordinate z=Y/X, a root corresponds to physical spinor ratio b/a=1/z. Zero is the south pole and infinity is the north pole. Missing leading powers give the multiplicity at infinity. Exact representatives use algebraic factorization or an exact Sylvester-matrix certificate for their multiplicities. Refined numerical representatives have separate interval root-count checks for the saved vectors; their AC equalities remain numerical. Sphere positions are numerical visualizations.

Full partitions generally describe more than one rotation orbit or SLOCC class. Invertible collective filters preserve multiplicities but generally change AC order.

## Reproducing the checks

The website needs only a browser. Rebuilding the data and symbolic verification requires Python 3 and SymPy (validated with SymPy 1.12). Numerical refinement additionally uses NumPy, SciPy and mpmath; the independent refinement audits use mpmath and the Python standard library. From this directory:

```bash
python3 research/verify_states.py
python3 research/acliege/verify_exact_states.py
python3 research/verify_multiplicities.py
python3 research/verify_moment_bounds.py
python3 research/partition-tetrahedra/verify.py
python3 research/n11-g2/verify.py
python3 research/n11-g2/verify_majorana.py
python3 research/n11-g2/audit_n11.py
python3 research/n11-g2/audit_n11_squarefree.py
python3 scripts/build_data.py
python3 research/ac-order-sparse.py
python3 research/ac-order-constructions.py
python3 research/ac-order-improvements-audit.py
python3 research/verify_partition_bounds.py
node scripts/build_symmetry_groups.cjs
python3 research/n12-four-double/verify-antipodal-obstruction.py
python3 research/n12-four-double/analytic-rotation-obstructions.py
python3 research/n12-four-double/analytic-d2-obstruction.py
python3 research/n12-four-double/verify-d2-independent.py
python3 research/n12-four-double/verify-boundary-and-rotation.py
```

The numerical audit commands also work from this catalogue directory:

```bash
python3 research/liege-refinement/audit/verify_refinement.py --original research/liege-library/observations.json --refined research/liege-refinement/refined.json
python3 research/liege-refinement/audit/verify_alternates.py --original research/liege-refinement/alternate-inputs.json --refined research/liege-refinement/alternate-refined.json
python3 research/liege-refinement/audit/verify_points.py --refined research/liege-refinement/refined.json --points research/liege-refinement/refined-points.json
```

To repeat the refinement, run `python3 research/liege-refinement/refine.py`; it reads frozen inputs beside the script. See the refinement README for alternate-vector and sphere-point regeneration commands. These programs preserve the originals and keep the target orders fixed.

The first-order existence proof is analytic; its claims are not inferred from numerical optimization. `research/partition-proofs.md` gives the argument. The partition verifier checks the combinatorial classification and bounds against the exact source records. `research/SMALL_N_PROOFS.md` and `research/ac-exclusions/index.html` describe the inherited AC exclusion certificates and their verifiers.

The scientific inputs live in `research/`. `scripts/build_data.py` and `scripts/partition_catalogue.py` generate `data/catalogue.json`, `data/catalogue.csv` and `data/catalogue.js`. The last format enables local-file browsing without fetch requests. Do not edit just one generated data format.

## GitHub Pages

Place this folder at the desired repository location and retain `.nojekyll`, `data/`, `research/`, `references/` and `vendor/`. No compilation is required. For companion navigation, keep `majorana-degeneracy-catalogue`, `anticoherent-state-catalogue`, `symmetric-sector-catalogue` and `ac-subspace-catalogue` as sibling folders.

This delivery is saved locally for future GitHub publication. KaTeX and its license are bundled under `vendor/katex/`; reference PDFs retain their original authorship and publication metadata.
