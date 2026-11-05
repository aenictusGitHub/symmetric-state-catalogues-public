# Pure anticoherent symmetric states

A static, offline catalogue indexed by qubit number **n** and anticoherence order **t**. Open `index.html` in a browser. KaTeX and the cited published source PDFs are included locally.

Each attached representative has a **View Majorana sphere** button. The interactive sphere supports dragging, keyboard rotation, multiplicity labels, symmetry axes, convex hulls, centre and barycentre markers, SVG export, and geometric point-group identification. Switching the representative selects its own constellation. Excluded or unresolved cells have no plotted witness. Coordinates are included when downloading a selected record and are also available in `data/majorana.json`. The representative panel shows its full point group and proper rotational subgroup without opening the sphere. AC representative choices and construction cards also include their group labels.

The sphere uses the same convention as the degeneracy catalogue: the homogeneous polynomial is `P(X,Y) = Σ_r (−1)^r √binom(n,r) c_r X^(n−r)Y^r`; a root `z = Y/X` represents the spinor ratio `b/a = 1/z`. Zero and infinity map to the south and north poles. Displayed coordinates and symmetry classifications are numerical; the existing exact AC certificates remain the evidence for anticoherence.

The main table covers every qubit number n=1–60 and orders t=1–14. The original Liège library supplies numerical witnesses through n=60, alongside analytical representatives and proved exclusions. Missing numerical solutions are left unresolved unless an independent proof excludes them.

Anticoherence means **at least** the selected order: the t-body reduced state is I/(t+1) on the symmetric subspace. Every multipole rank 1,…,t must vanish. The condition differs from maximal mixing on the full 2^t-dimensional qubit space and from cancelling only even multipoles.

Where an exact-order representative is recorded, it is the default for the selected cell. Higher-order states remain available as alternatives. Analytical exact-order witnesses are checked by two exact reduced-state calculations, including a nonzero deviation at the next order. Numerical library states are labelled separately and display reduced-state residuals rather than an exact existence claim. Read `research/exact-orders/index.html` for the construction and saved certificates.

A general rational construction and proof are given in [the exact-order research note](research/exact-orders/general/index.html). It guarantees exactly order t for every N ≥ 4t(t+1)^4; the bound is sufficient, not optimal. A separate three-component family gives exactly order two for every N ≥ 7.

The highest primary Julia/website order at each size is labelled **Conj. max** unless exact attainment and a global upper bound establish **Proved max**. These labels concern the maximum order at the qubit number; each representative retains its own exact or numerical evidence. Lower-order and superseded library states are not marked as current global maxima.

## Evidence

The original core at n=2–12, t=1–6 comprises 26 exact existence results and 30 proved exclusions. The expanded table contains exact witnesses, numerical witnesses, proved exclusions and unresolved cells. A failed numerical search is never treated as a proof of impossibility.

The analytical representatives have their normalization and every reduced matrix through the claimed order checked symbolically. The dataset includes multiple representatives at some sizes and independently derived radical states at n=11, 14 and 17. A witness proves existence, not the largest achievable order.

For n=2,…,12, the highest orders witnessed here are 1,1,2,1,3,2,3,2,3,3,5. Maximality is established for every n=2–12. In particular, n=11 has maximal order 3, and n=12 is the first qubit number admitting pure fourth-order anticoherence.

The previously reported exclusions (9,3) and (10,4) now have independent exact proofs: a four-copy positivity certificate and a three-copy polynomial identity, respectively. The n=11,t=4 exclusion is also now proved by both a 26-term polynomial identity and an invertible 27-by-27 matrix of necessary three-copy equations (determinant 3 modulo 37). Read `research/reported-exclusions/index.html` for the proofs, LaTeX formulas, and reproducible certificates.

Two-copy exact identities exclude (8,4) and (10,5). Schmidt-rank, odd-complement, and tensor-Gram arguments establish the other exclusions. Read `research/SMALL_N_PROOFS.md` for the derivations.

The expansion imports the audited library representatives through n=60, with independently checked high-precision partial traces. The current dataset and [import and audit notes](research/library42/index.html) give representative counts, residuals and rejected source candidates, including the starred N=29,t=6 candidate and the mislabelled acs.48.8.txt file. The research folder retains its original name to keep existing links valid.

## Reproduce

With Python 3 and SymPy, from this folder:

```sh
python3 research/verify_states.py
python3 research/verify_small_n_audit.py
python3 research/reported-exclusions/n9-t3/verify.py
python3 research/reported-exclusions/n9-t3/independent_verify.py
python3 research/reported-exclusions/n10-t4/verify.py
python3 research/reported-exclusions/n11-t4/verify.py
python3 research/reported-exclusions/n11-t4/moment/verify.py
python3 research/reported-exclusions/n11-t4/moment/independent_verify.py
python3 research/exact-orders/build.py  # requires NumPy and SymPy
python3 research/exact-orders/verify.py
python3 research/exact-orders/build_page.py
python3 research/library42/build.py
python3 research/library42/verify_exact.py
python3 research/library42/verify_numeric.py
python3 research/library42/build_page.py
python3 scripts/build_data.py
```

The first verifies all representative states in the supplied file; pass `data/catalogue.json` to include every representative in the rebuilt dataset. The second independently verifies the core witnesses and both exact rational exclusion identities. The remaining proof verifiers check all new higher-copy exclusions using exact arithmetic. The last rebuilds browser data and CSV from the source records.

After changing representative data, regenerate the Majorana coordinates with `python3 scripts/build_majorana.py` (SymPy and NumPy). Independent reconstruction and viewer checks are included:

```sh
python3 research/majorana-viewer/verify_coordinates.py
node research/majorana-viewer/verify_viewer.cjs
node research/majorana-viewer/verify_app.cjs
```

The coordinate verifier uses mpmath; the JavaScript checks use the locally bundled KaTeX. Saved reports are in `research/majorana-viewer/`.

`data/catalogue.json` is the complete distributable dataset; `data/catalogue.js` is its offline browser version. `research/states.json` and `research/restrictions.json` preserve witness definitions and source audits. Downloads contain complex amplitudes as [real, imaginary] pairs and exact symbolic coefficients separately.

## Source corrections and scope

PRResearch 4, 013178 (2022), Table IV, has two inconsistent rows: its nine-qubit state is only second-order anticoherent, and its printed ten-qubit coefficients are not normalized. Exact reduced-state checks and the independently checked coefficient correction are recorded in `research/table-iv-audit.json`. These are catalogue corrections, not claimed published errata.

The 2017 numerical survey's missing points indicate failed convergence. The cyclic-construction thresholds in the 2015 paper apply to a restricted ansatz and are not universal thresholds. The site preserves these distinctions and links directly to the relevant PDF pages.

The companion sector-length catalogue is linked as a sibling folder. Keep both catalogue folders beside each other when publishing them together on GitHub Pages. This site has no server-side dependencies and no runtime network requests.

The supplied Julia formulas add exact representatives at n=13,15,16,18,20 and simpler alternatives at smaller sizes. See `research/acliege/pure-state-audit.json` and its portable verifier. The companion [AC-subspace catalogue](../ac-subspace-catalogue/index.html) includes full-span certificates and new radical reconstructions.

## Cross-check with the AC-subspace catalogue

The two catalogues share compatible existence and exclusion statements. The cross-check imports the already certified fourth-order radical states at n=14 and n=17; `research/cross-catalogue/import_exact.py` independently verifies every reduced matrix through order five. Both imported representatives have exactly order four. The later three-copy exclusion proves this is globally maximal at n=14; maximality at n=17 remains unresolved. Their Majorana coordinates use the source's exact square-free factorization and are checked by reconstructing the Dicke vector from the imported points.

The same update applies the general rectangular-tensor bounds and their exact equality obstructions. Former unknown cells (15,6), (16,6), (16,7), (18,7) and (20,7) are now proved impossible. In particular, the sixth-order sixteen-qubit exclusion proves that the existing fifth-order sixteen-qubit state has globally maximum AC order. Full proofs and independent certificates are linked to the sibling AC-subspace catalogue.

The local Liège library contains a thirteen-qubit third-order Dicke file, with a nonzero fourth-order residual. No thirteen-qubit fourth-order Dicke file was found. A missing file is not a nonexistence proof. The notation `acs.N.t.g.txt` uses its last index for Majorana-point multiplicity, not AC order; it must not be read as a t=g witness.

The thirteen-qubit fourth-order case is now separately proved impossible by the exact four-copy certificate in `../ac-subspace-catalogue/research/n13-t4/index.html`. Its 320 independent necessary equations on 331 real invariant coordinates leave an eleven-dimensional kernel, and the physical trace functional vanishes identically on that kernel. A normalized four-copy moment has trace one, giving a contradiction. Both former unknown cells (13,4) and (13,5) are now proved impossible, and the existing thirteen-qubit third-order representative is globally maximal.

The exact three-copy positivity certificate in `../ac-subspace-catalogue/research/n14-t5/index.html` now proves that no fifth-order fourteen-qubit pure state exists. The known fourth-order radical representative is therefore globally maximal. The catalogue and independent cross-check include this new exclusion.

## Fifteen-qubit maximum order

The maximum pure-state AC order is **three**. Fourth and higher orders are excluded by an exact four-copy partial-transpose positivity certificate, independently verified in rational arithmetic. See [the full proof](../ac-subspace-catalogue/research/n15-t4/index.html). The existing exact Julia-source third-order state attains the maximum.

## Seventeen-qubit maximum order

The maximum pure-state AC order at N=17 is **four**. The existing analytical state attains it; a [three-copy exact certificate](../ac-subspace-catalogue/research/n17-t5/index.html) excludes every fifth-order state. Sixth order also has an independent [rank-one-complement proof](../ac-subspace-catalogue/research/n17-t6/index.html), part of a general infinite family of exclusions. All higher orders are excluded by nesting.

## Eighteen-qubit maximum order

The maximum pure-state AC order at N=18 is **five**. The existing exact four-component Julia-source state attains it. An independent exact [three-copy trace–probability certificate](../ac-subspace-catalogue/research/n18-t6/index.html) excludes every sixth-order state and all higher orders. The [fifth-order subspace maximum](../ac-subspace-catalogue/research/n18-t5/index.html) is one.

## Nineteen-qubit maximum order

The maximum pure-state AC order at N=19 is **four**. A new [five-component radical state](../ac-subspace-catalogue/research/n19-t4/index.html) upgrades the previous numerical existence report. An independently verified [four-copy matrix-positivity certificate](../ac-subspace-catalogue/research/n19-t5/index.html) excludes every fifth-order state, and a separate [three-copy probability identity](../ac-subspace-catalogue/research/n19-t6/index.html) excludes sixth order. Both cover arbitrary complex amplitudes.

The eighteen-qubit AC5 representative also spans an optimal fourth-order AC subspace: [K(18,4)=1](../ac-subspace-catalogue/research/n18-t4/index.html), by an exact auxiliary-spin-three inertia proof. This is a subspace-dimension result; its maximum pure AC order remains five.


## Twenty-qubit maximum order

The [exact three-copy matrix-positivity certificate](../ac-subspace-catalogue/research/n20-t6/index.html) excludes every pure sixth-order AC state on twenty qubits, including arbitrary complex amplitudes. Its independent verifier reconstructs 95 selected physical equations and proves all 61 rational factorization pivots positive. Consequently all orders t≥6 are impossible. The existing exact Julia-source AC5 state attains globally maximum pure-state order five. The N=19 AC4 representative also links to the separate four-copy proof that its one-dimensional span has maximum fourth-order subspace dimension.
