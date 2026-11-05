# No pure symmetric nine-qubit state is 3-anticoherent

This is an exact finite-dimensional certificate. Its arithmetic is integer/rational; numerical searches are not part of the proof. The certificate applies to arbitrary complex amplitudes, not just real states or states with a chosen symmetry.

Let V be the spin-9/2 space, of dimension 10. If a pure vector psi is 3-anticoherent, every spherical-tensor expectation of ranks 1, 2 and 3 vanishes. This property is invariant under SU(2) rotations.

## 1. A hypothetical state gives a positive four-copy moment

Twirl the fourth tensor power of its projector:

X = integral over SU(2) of (U |psi><psi| U†) tensor-powered four times.

Then X is nonzero, positive semidefinite, SU(2)-invariant, and supported on Sym^4(V), of dimension 715. For every rank-ell tensor, ell=1,2,3, it obeys

Tr_1[(T_(ell,0) tensor I^3) X] = 0.

Indeed, before integration this partial contraction is the zero tensor expectation multiplied by the third tensor power of the rotated projector. We only use q=0 necessary equations.

The spin decomposition is

Sym^4(V_(9/2)) = 2 V_0 + 3 V_2 + 2 V_3 + 4 V_4 + 2 V_5 + 5 V_6 + 3 V_7 + 4 V_8 + 3 V_9 + 4 V_10 + 2 V_11 + 3 V_12 + V_13 + 2 V_14 + V_15 + V_16 + V_18.

For each spin J, write its invariant block as a Hermitian multiplicity matrix B_J. Each B_J is positive semidefinite. A real Clebsch–Gordan/binomial basis is used. Entrywise real parts of these Hermitian matrices remain real positive semidefinite, remain nonzero collectively, and preserve all necessary equations because the equations have real coefficients. Thus it suffices to rule out real symmetric positive blocks. Their total number of independent entries is 88.

## 2. Rational convention fixes every matrix entry

For an occupation vector o of four factors, use E_o, the *averaged* symmetrized tensor product of b_i = sqrt(binomial(9,i)) |D_9^(i)>. Its squared norm is

  product_i binomial(9,i)^(o_i) / multinomial(4;o).

The collective raising coefficient o -> o-e_i+e_(i+1) is o_i(i+1), and the lowering coefficient is o_i(10-i). Hence highest-weight vectors and their descendants can be represented rationally. The complete highest-weight matrices H_J are stored in kernel.json and checked by the verifier.

Put U_(J,s) = (J_-)^s H_J/s!, for 0 <= s <= 2J. The invariant contribution is

  d_J sum_s U_(J,s) B_J U_(J,s)^T / binomial(2J,s),

Here the products U B U^T denote sums of dyadic outer products of the actual descendant vectors, with their physical inner product; they are not ordinary operator-coordinate products in the nonorthonormal occupation basis. The number d_J is the positive integer denominatorScale in kernel.json. These positive factors only fix the scale of B_J and do not affect positivity. The ladder relations are verified exactly. In the stored JSON, J2 means 2J and multiplicity indices start at zero.

The three diagonal tensors may be scaled to the integer polynomials in z=2i-9:

  f_1(i)=z,   f_2(i)=3 z^2-99,   f_3(i)=z(5 z^2-293).

Taking the first factor against basis state |i> sends E_o to sqrt(binomial(9,i)) o_i E_(o-e_i)/4. Thus every contraction equation has rational coefficients. A common factor 1/16 is omitted in the equations.

## 3. An exact positive diagonal identity forces zero

Let b be the vector of the 88 real block entries, in the order variables in kernel.json. Eighty-five specified necessary contraction equations form A b=0. The portable verifier rebuilds these equations directly from the representation above. It checks rank(A)=85 and supplies an exact 88-by-3 kernel basis K. Consequently every feasible b equals K y for some real y.

Define the four diagonal entries, in the fixed highest-weight bases:

  b0 = (B_0)_(0,0),
  b1 = (B_2)_(0,0),
  b2 = (B_2)_(1,1),
  b3 = (B_3)_(1,1).

Every b in ker(A) obeys the exact identity

  2003631119000 b0 + 171013438791 b1
  + 5530657636 b2 + 12265113288000 b3 = 0.

All four coefficients are strictly positive. Since the B_J are positive semidefinite, all four diagonal entries are nonnegative. They must therefore all vanish. The verifier checks that the first three corresponding rows of K are linearly independent. Hence b0=b1=b2=0 forces y=0, and all B_J vanish. This contradicts the nonzero four-copy moment X.

Therefore no pure symmetric nine-qubit state is 3-anticoherent. Higher orders are excluded as well, since higher anticoherence implies order 3.

## Reproduction

Run `python3 verify.py` in this folder (or invoke its full path). It requires SymPy and reads only its neighboring kernel.json and certificate.json. It performs exact reconstruction and rank/identity checks; NumPy, SciPy, SDP software and floating-point arithmetic are not required. The independent verifier uses a separate symmetrized-word derivation of the contraction coefficients.

Required proof files: verify.py, kernel.json, certificate.json. proof.md and formula.tex provide the mathematical explanation. verification.json records a successful run but is not trusted as proof input.
