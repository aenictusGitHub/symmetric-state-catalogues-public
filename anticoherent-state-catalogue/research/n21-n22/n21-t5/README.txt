N=21: maximum pure-state anticoherence order is four

This bundle proves nonexistence at order five using a six-copy separating
functional, and checks an explicit state of exact order four.

Requirements: Python 3, NumPy and python-flint (Arb). SciPy is not needed.
Run from this directory:

    python3 verify_certificate.py

The verifier recomputes every representation and positivity bound using
exact integer factorial sums and outward-rounded Arb arithmetic at 128 bits.
All values stored in the NPZ input matrices are interpreted as exact binary
rationals. They are certificate inputs, not trusted numerical conclusions.

A successful run writes independent-verification.json and
bootstrap-verification.json, including exact rational bounds and SHA-256
hashes. The expected run time is approximately two minutes on the machine
used to prepare this certificate.

proof.tex gives the mathematical argument, normalization and perturbation
bounds. No discovery scripts, numerical optimization or large moment arrays
are required.
