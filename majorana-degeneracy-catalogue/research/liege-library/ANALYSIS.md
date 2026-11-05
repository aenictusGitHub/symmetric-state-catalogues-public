# Numerical audit of selected Liège-library counterparts

This analysis contains numerical observations, not new exact anticoherence or degeneracy proofs. Only filenames indexed by the two primary library pages were selected. The website host was unreachable, so local counterpart files were used with recorded paths and hashes. Identity with the website's current bytes remains unverified.

The final audit covers 144 of 145 indexed records: 99 general states and 45 of 46 degenerate-library states. The coordinate counterpart acs.12.2.3.txt was unavailable. It is recorded as missing; no numerical state or partition is inferred for it.

All 144 available records satisfy the chosen numerical consistency criterion: squared Frobenius distance of rho_t from I_(t+1)/(t+1) at most 10^(-20). The largest observed value is 8.21464475140586e-27, for acs.52.9. This check supports the source's claimed order numerically; it is not a proof of exact order or maximality.

For the 99 records with both coordinate and Dicke files, 83 pairs match after a global phase and explicit basis-convention transformations (reversal, conjugation, alternating signs). Their largest infidelity is 8.89932204914457e-52. Those records use normalized original Dicke coefficients after the matching transformation. The remaining 16 records use coefficients reconstructed directly from their own coordinate files, ensuring the displayed state and plotted points refer to the same numerical representative. Failure of the tested matching transformations does not prove physical inequivalence: arbitrary spatial rotations were not optimized, and different local versions are also possible.

The 16 unmatched pairs are acs.2.1, acs.3.1, acs.4.2, acs.5.1, acs.6.3, acs.7.2, acs.8.3, acs.9.2, acs.10.3, acs.11.3, acs.12.5, acs.13.3, acs.15.3, acs.16.5, acs.18.5 and acs.20.5.

Every inferred full multiplicity partition is unchanged when the coordinate-clustering tolerance varies over 10^(-12),10^(-10),10^(-8),10^(-7),10^(-6),10^(-4). All repeated clusters have zero numerical diameter, while the smallest separation between distinct clusters is 0.08317724931665527, for acs.11.3.2. Thus there is a large observed separation between coincident rows and distinct points, but the observations remain labelled numerical.

A useful new core observation is acs.12.3.2: its inferred full partition is(2,1^10), rather than the(2^6) partition of the exact double-octahedral representative already in the catalogue. Its repeated south-pole rows are identical, its nearest distinct-point separation is 0.3917648403039392, and its reconstructed rho3 squared Frobenius error is 9.85836956554273e-41. This is a source-reported numerical candidate in a new partition cell, not an exact lower-bound upgrade.

## Method and reproduction

`analyze_records.py` reads only manifest-selected files. It accepts both the Mathematica complex-number format with precision marks and the real/imaginary two-column format present in the selected counterpart files. No unindexed local states are scanned or imported.

Each Cartesian point is first normalized to the unit sphere; the original radius error is recorded. Its physical spinor is formed using a north or south coordinate chart chosen to avoid a small denominator. The product polynomial of these spinors is accumulated at 90 decimal digits. If its degree-r coefficient is p_r, the normalized Dicke coefficient is proportional to p_r/sqrt(binomial(n,r)). The resulting reduced matrix is computed directly with the exact combinatorial factors at the same numerical precision.

`observations.json` contains the full per-record evidence: coefficients, plotted points, high-precision strings, norm errors, reduced-matrix residuals, all clustering tolerances, coordinate/Dicke comparisons and source provenance. The normalization and clustering steps do not erase the source caveats.

Run `python3 analyze_records.py` with mpmath installed. Both script paths and raw-file paths can be relative to this directory. The batch takes approximately 12–15 seconds on the development machine. `validate_analysis.py` separately checks the reconstruction on locally constructed tetrahedral and octahedral constellations, their repeated versions, and a 100-point regular polygon; these validation examples are explicitly not labelled imported source data.
