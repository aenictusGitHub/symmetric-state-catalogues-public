AC maximum-order records compared with spherical designs on S^2

The reference is the Hardin–Sloane table of putatively optimal exact-N,
equal-weight spherical-design strengths in ambient dimension three:
  https://neilsloane.com/sphdesigns/index.html
The associated paper is Discrete Comput. Geom. 15 (1996), 429–441:
  https://doi.org/10.1007/BF02711518
The table says last modified 2 October 2002; retrieved 6 November 2026.
Cohn's maintained mirror also calls these strengths putatively optimal:
  https://cohn.mit.edu/sloane/

This compares exactly N distinct equal-weight points on S^2 with the
reported maximum AC order for N symmetric qubits. It does not measure the
design strength of any particular state's Majorana constellation. Existence
at fewer points is not treated as existence at exactly N points.

Current recorded-order gaps:
  N=29: reported AC maximum 5, design benchmark 6.
  N=87: reported AC maximum 11, design benchmark 12.
  N=94: reported AC maximum 12, design benchmark 13.
The other 97 recorded orders agree, including the vacuous order-zero case
at N=1. This is a comparison of records, not a proof about true maxima.

AC maximum orders at all three gaps remain conjectured. At N=29 an exact
fifth-order witness and a seventh-order exclusion give rigorous bounds 5–6.
At N=87 and N=94, the reported AC orders are numerical witnesses. The design
source marks all three comparison entries V3, meaning numerical coordinates
with source-reported discrepancy at most 1e-26, rather than a V1 algebraic or
V2 interval existence proof. These evidence codes concern existence of the
construction, not optimality. No design-based AC upper bound is introduced.

Files
  benchmark.json / benchmark.csv: factual N/order/evidence fields, with source
    attribution and a hash of the retrieved source HTML (the HTML is not hosted).
  comparison.csv: all 100 reported-order comparisons and rigorous AC bounds.
  comparison-audit.json: independent checks and the three reported-order gaps.
  import_benchmark.py: standard-library importer for a saved copy of the source.
  verify_comparison.py: independent audit using the published catalogue data.

Rebuild the factual benchmark after saving the source HTML:
  python3 import_benchmark.py --html saved-spherical-designs-source.html
Audit the comparison:
  python3 verify_comparison.py
An optional --baseline catalogue.json asserts that all existing scientific
states, existence entries, proof records and proposed AC maxima are unchanged.
