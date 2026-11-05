# Refined anticoherent-state components

This archive contains 144 primary coefficient files in `coefficients/` and 16 separate original-Dicke alternatives in `alternates/`.

Each file is plain text with n+1 rows. Row r (starting with r=0) is the coefficient of the normalized Dicke state |D_n^(r)>. The first column is its real part and the second its imaginary part, separated by a tab. Nonzero components are exported to 85 significant digits. There are no header rows.

The `acs.*` source identifiers and originally reported AC orders are retained. All primary vectors have a numerically verified Hilbert-Schmidt residual below 1e-60 at that reported order; the largest is 5.872e-66. The 16 alternate Dicke vectors also meet that target. These are numerical refinements, not new exact AC proofs or searches for a higher order.

Original files are unchanged. The separate research package includes original/refined JSON arrays, refined sphere points, a download manifest with hashes, reproduction scripts and independent audits. Independent root-count checks establish the full partitions of the saved primary numerical vectors; alternates are not assigned a partition or constellation.

Nine primary states already met the target and were unchanged. Five singular second-order cases use their forced two-term Dicke family, removing tiny forbidden amplitudes while retaining the original nonzero phases. The remaining corrections use high-precision residuals and an analytic Jacobian.
