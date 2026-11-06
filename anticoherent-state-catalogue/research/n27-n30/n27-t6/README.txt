Pure order-6 anticoherence on 27 symmetric qubits: six-copy certificate

This bundle verifies a global exclusion for arbitrary complex pure states.
The exact order-5 radical witness in the surrounding n27-n30 research folder
then makes the maximum AC order equal to 5.

The mathematical theorem depends on a COMPLETE successful verification run.
A numerical separator, solver log, or partly completed check is insufficient.

Reproduction
  Requires Python 3, NumPy, and python-flint (Arb). SciPy is not required.
  Keep every file listed in MANIFEST.txt together, then run:
      python3 verify_certificate.py

The verifier reconstructs the physical symmetric-power representations,
all used AC contraction coefficients, the first-leg partial-time-reversal
adjoint, and every source positivity bound using exact integer factorial
sums and outward-rounded 128-bit Arb arithmetic. All archive floats are
interpreted as exact binary rationals. The PPT duals are positive by exact
Gram factors. No optimizer output or floating residual is accepted as a
proof bound. The output sectors used are doubled spins 0,2,...,20 for the
1|5 split; conditional Kramers/Breuer-Hall cuts are not part of this proof.

Outputs
  bootstrap-verification.json: validated recursive CFP construction.
  independent-verification.json: written only after all 81 source sectors
    pass strict positivity, with rational bounds and all input SHA256 hashes.

proof.tex supplies definitions, normalizations, the separating identity,
the representation bootstrap, and the complete positivity argument.
The large CFP and dual archives are split into files below GitHub's limit.
Discovery scripts are not needed for independent verification.
