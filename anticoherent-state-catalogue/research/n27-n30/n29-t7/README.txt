Validated pure-state AC exclusion: N=29, order=7.

Requires Python3, NumPy and python-flint (Arb); SciPy is not required.
Run python3 verify_certificate.py. All NPZ floats encode exact binary rational
inputs. Every representation, contraction and positivity bound is recomputed
with exact integer factorial sums and outward-rounded128-bit Arb arithmetic.
The run writes independent-verification.json and bootstrap-verification.json.
proof.tex gives the mathematical argument.

This certificate also uses selected positive2|2 partial-time-reversal blocks.
Their dual matrices are exact dyadic Gram factors, hence positive by construction.
The third-copy CFP perturbation in pair-pair recoupling is explicitly included.
