# Refinement of the reported anticoherent states

The 144 primary numerical representatives were refined at their original reported orders. All now have a Hilbert–Schmidt residual below 10^(-60), verified after reading back the exported coefficients. The largest final residual is 5.871093505912389e-66. Nine exact balanced Dicke representatives already met the target and were unchanged; the other 135 were improved. No higher-order search was performed and no numerical representative was promoted to an exact proof.

The 16 original Dicke counterparts that did not match their coordinate representatives under the tested convention changes were refined separately. All 16 also meet the target, with largest residual 1.333085437680919e-70. They remain download-only alternates: they have no assigned full partition or associated sphere plot in this package.

Every original file is unchanged. The frozen inputs, original coefficient arrays, and before/after residuals are retained.

## Constraints and precision

For the reported order t, form the symmetric Schmidt matrix

  A_(a,q) = c_(a+q) sqrt[binomial(t,a) binomial(n-t,q) / binomial(n,a+q)].

The equations are AA†=I_(t+1)/(t+1), with each off-diagonal real and imaginary residual weighted by sqrt(2). The Euclidean norm of this real residual vector is therefore exactly the Hilbert–Schmidt norm of rho_t-I/(t+1). These unnormalized equations include the state-normalization constraint through their trace.

Each ordinary step uses the analytic real Jacobian, a double-precision SVD approximate right inverse, and 100-digit residual evaluation and corrections. The output has 85 significant digits per real and imaginary component. The reported final residual is recalculated from those serialized values. Printed precision and a small numerical residual are not assertions of an exact algebraic solution.

For repeated south-pole points, c_0 through c_(g-1) are held identically zero. For repeated north-pole points, c_(n-g+1) through c_n are held identically zero. These structural constraints prevent repeated poles from splitting. Simple points may move.

## Five singular second-order families

For n=3g+1 and a repeated pole of multiplicity g, second-order moment saturation forces a two-term Dicke family. The affected source records are n=7,10,13,16,19 with g=2,3,4,5,6.

For a south pole, the only nonzero amplitudes are c_g and c_n, with squared magnitudes n/(4g+2) and (g+1)/(4g+2). The north-pole case is reversed. The correction removes tiny forbidden interior amplitudes, sets the two required magnitudes, and retains their original phases. This is the nearest vector in that forced family in Euclidean coefficient distance. It should not be described merely as changing the number of printed digits.

The largest component change among all primary refinements is 1.135802503137318e-7, at n=19. Most other changes are much smaller. The family corrections preserve the original full partition.

## Refined sphere points and independent checks

The refined sphere points are matched to the refined coefficients, rather than reusing the original plot without a check. Each simple Majorana root was corrected from its source point by 110-digit Newton iteration in a well-conditioned direct or reciprocal chart. Repeated poles remain fixed. All 144 root passes succeeded.

The largest point displacement is 1.801215599248130e-7. The minimum distance between distinct refined points remains 0.08317724928979622. The original point order is retained in the high-precision arrays.

An independent implementation checks norms, reported-order reduced matrices, component changes, exact pole zeros and the serialized precision. It also establishes one simple root in each of 5,437 disjoint neighborhoods, together with 238 fixed pole roots. This confirms that the full partitions are preserved for the exported primary vectors. These root-count checks do not turn the approximate AC equalities into exact AC proofs.

## Files and reproduction

- `original-states.json`: frozen high-precision inputs for the 144 primary representatives.
- `refined.json`: original and refined coefficients, residual history, normalization errors, structural zeros and method details.
- `refined-points.json`: refined sphere points, point movements and numerical matching data, separate from the coefficient audit file.
- `alternate-inputs.json` and `alternate-refined.json`: the 16 unmatched original Dicke counterparts, with no partition claim.
- `refine.py`: the refinement procedure; requires NumPy, SciPy and mpmath.
- `refine_points.py`: recomputes sphere points from the exported primary coefficients.
- `audit/`: independent verification scripts and reports.

Run `python3 refine.py` to repeat the primary refinement. Run `python3 refine.py --input alternate-inputs.json --output alternate-refined.json` for the alternates, then `python3 refine_points.py` for primary sphere points. Paths default to files beside the scripts. The two optional prepare scripts recover the frozen inputs from the neighboring `liege-library` package; they do not search for additional states.

## Downloads and independent verification

`coefficients/` contains 144 primary files named `{id}.refined.txt`; `alternates/` contains 16 separate alternate files. Each has n+1 rows with tab-separated real and imaginary Dicke coefficients, in increasing excitation order r=0,...,n. `download-manifest.json` records every path, row count and SHA-256 hash.

`refined-components.zip` contains exactly these 160 coefficient files and a format README. It does not include itself, logs or caches. The complete research files remain available alongside the ZIP.

From this research directory, run:

```bash
python3 audit/verify_refinement.py --original ../liege-library/observations.json --refined refined.json
python3 audit/verify_alternates.py --original alternate-inputs.json --refined alternate-refined.json
python3 audit/verify_points.py --refined refined.json --points refined-points.json
```

These independent checks use mpmath and the Python standard library. The full primary audit includes interval root-count checks; it takes longer than the coefficient-only and sphere-point consistency checks.
