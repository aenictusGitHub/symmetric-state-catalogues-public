> Current status: the separate auxiliary-spin inertia theorem proves k_max(13,2)=2. This file records an earlier relaxation that did not decide existence.

# Three-copy projector relaxation for N = 13, t = 2, k = 3

This investigation **does not settle existence**. It strengthens the earlier single exterior-power relaxation by retaining all three permutation sectors of the moment of a hypothetical rank-three projector. The saved candidate satisfies the numerical equations and positivity conditions to floating-point accuracy. It is not an actual subspace and is not an exact feasibility certificate.

The calculation has 168 real invariant variables, 70 independent homogeneous equations, and 352 positive-semidefinite blocks. The recorded solve returned `optimal_inaccurate`; direct numerical checking found equality residual below 2 × 10⁻¹⁷, trace within 2 × 10⁻¹² of one, and minimum eigenvalue approximately −1.1 × 10⁻¹⁴. These values describe numerical consistency, not a proof.

`proof.tex` explains the necessary conditions and their operator ordering. `candidate.json` stores the candidate in a fixed Clebsch–Gordan basis, independently of the arbitrary bases used to parameterize permutation sectors. The construction keeps all symmetric, antisymmetric, and mixed-symmetry sectors; it is distinct from the earlier exterior-power-only test.

## Reproduce without an optimizer

Install NumPy and SciPy, then run:

```sh
python generate.py 13 3
python verify_candidate.py
python audit_positive_controls.py
```

Generated maps and fresh reports go into `generated/`. Set the environment variable `AC_P3_OUTPUT` to choose another output folder. Large generated matrices are intentionally omitted from the distributed catalogue; the generator reconstructs them.

The controls populate the moment from actual known codes at N = 7, k = 2 and N = 14, k = 3, and compare the maps with direct tensor products. The stored N = 14 basis is the numerical evaluation of the catalogue's exact algebraic witness. Its exact construction and certification are in [the N = 14 research folder](../../algebraic-k3/n14/).

## Optional new numerical solve

Install CVXPY and CLARABEL in addition to NumPy and SciPy, then run:

```sh
python solve.py 13 3
```

The solver has a 90-second optimization limit; preparing the semidefinite problem also takes time. Running the verifier against the saved `candidate.json` does not require CVXPY. Solver feasibility, numerical tolerances, and failed searches are never used here as an exact existence or exclusion proof.
