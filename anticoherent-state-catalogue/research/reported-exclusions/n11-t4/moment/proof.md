# No pure symmetric eleven-qubit state is 4-anticoherent

This exact certificate applies to arbitrary complex pure states in the symmetric eleven-qubit space. Numerical optimization is not used in the proof.

Let V be the spin-11/2 representation. Suppose a normalized vector psi in V were 4-anticoherent, so its spherical-tensor expectations of ranks 1 through 4 all vanish. Form its third-copy SU(2) moment

  X = integral over SU(2) of (U |psi><psi| U†) tensor-powered three times.

This is a nonzero positive invariant operator on Sym^3(V), with trace one. Anticoherence is preserved by rotations, so

  Tr_1[(T_(ell,0) tensor I^2) X] = 0, for ell=1,2,3,4.

Each integrand vanishes because the first factor is a zero tensor expectation. Only these necessary q=0 contraction equations are used.

## Invariant blocks and complex states

The exact decomposition is

Sym^3(V_(11/2)) = V_(3/2) + V_(5/2) + V_(7/2) + 2V_(9/2) + 2V_(11/2) + 2V_(13/2) + 2V_(15/2) + 2V_(17/2) + V_(19/2) + 2V_(21/2) + V_(23/2) + V_(25/2) + V_(27/2) + V_(29/2) + V_(33/2).

Its dimension is 364. Choose real highest-weight bases. The invariant operator is represented by Hermitian positive multiplicity matrices B_J. Entrywise real parts of those matrices remain real symmetric positive matrices and preserve all contraction equations, whose coefficients are real. They cannot all vanish for a nonzero positive operator. It therefore suffices to rule out nonzero real symmetric invariant solutions. Their total number of independent entries is 27.

## Rational basis and contraction equations

Use b_i = sqrt(binomial(11,i)) |D_11^(i)> and the averaged symmetric occupation basis E_o for three factors. Its squared norm is

  product_i binomial(11,i)^(o_i) / multinomial(3;o).

The collective raising and lowering coefficients are integer: o_i(i+1) for i to i+1, and o_i(12-i) for i to i-1. Rational highest-weight matrices H_J are stored in certificate.json. The verifier checks the highest-weight equations, independence, completeness, and the exact ladder relations.

With U_(J,s) = (J_-)^s H_J/s!, each invariant block has the form

  d_J sum_s U_(J,s) B_J U_(J,s)^T / binomial(2J,s).

The positive integer d_J is denominatorScale in the data. It only fixes the scale of B_J. All multiplicity indices are zero-based; J2 in the data means 2J.

For z=2i-11, the four diagonal tensor components can be scaled to

  f1=z,
  f2=3z^2-143,
  f3=z(5z^2-425),
  f4=35z^4-4190z^2+57915.

These are the degree-one through degree-four diagonal irreducible tensor polynomials. The verifier additionally checks their pairwise trace orthogonality, including the rank-zero constant.

Contracting the first factor against |i> sends E_o to sqrt(binomial(11,i)) o_i E_(o-e_i)/3. This gives rational coefficients for every necessary contraction equation. The common factor 1/9 is omitted. Twenty-seven specified equations are stored by rank, output weight and output matrix-entry labels; the verifier reconstructs them directly from these definitions rather than trusting the supplied matrix.

## Exact contradiction

Writing the 27 real block entries as a vector b, the chosen contraction equations give A b=0 for a 27-by-27 integer matrix A. Its determinant obeys

  det(A) = 3 modulo 37.

In particular, A is invertible over the rational and real numbers. Hence b=0, contrary to the nonzero invariant moment. No pure symmetric eleven-qubit state is 4-anticoherent, and higher anticoherence orders are also impossible.

## Reproduction

Run `python3 verify.py` from any directory, using its full path if needed. It reads the adjacent certificate.json and uses Python's standard library plus SymPy. All checks are exact integer/rational calculations; runtime is about 0.2 seconds. The independent verifier, `independent_verify.py`, reconstructs the same equations separately using symmetrized-word orbit multiplicities and imports no discovery or primary-verifier code.

Required proof files: verify.py and certificate.json. proof.md and formula.tex give the mathematical explanation. The two verification JSON files record successful runs but are not inputs trusted by either verifier.
