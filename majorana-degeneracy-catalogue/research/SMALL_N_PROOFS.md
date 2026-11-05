# Small-n pure anticoherence: independent audit

Full t-anticoherence means all ranks1 through t vanish, not merely the even ranks.

## Schmidt-rank bound

A pure symmetric n-qubit t-anticoherent state requires t≤floor(n/2).

The t-body reduced state is I_sym^(t)/(t+1), of rank t+1. Its complementary symmetric subsystem has dimension n−t+1. Equality of the nonzero Schmidt spectra requires t+1≤n−t+1. Higher orders are excluded by nesting.

## Odd-size complement obstruction

No pure symmetric n=2t+1 qubit state is t-anticoherent when t≥2.

For a hypothetical state, rho_t=I_sym^(t)/(t+1). The complementary rho_(t+1) has rank t+1 with equal nonzero eigenvalues in dimension t+2, hence rho_(t+1)=(I_sym^(t+1)−|phi><phi|)/(t+1). Since Tr_1 I_sym^(t+1)=(t+2)/(t+1) I_sym^(t), tracing one qubit forces Tr_1|phi><phi|=I_sym^(t)/(t+1). This reduced state has rank t+1>2, impossible for a pure (t+1)-qubit phi across the t-versus-one partition.

## Orthogonal tensor-image dimension bound

Pure t-anticoherence with t≥2q requires (q+1)^2≤n+1. In particular, n=12,t=6 is impossible.

Take all orthonormal irreducible tensors T_LM with 0≤L≤q. Products T_LM†T_L′M′ have ranks at most 2q≤t. Anticoherence makes their expectations equal to those in I/(n+1), so the (q+1)^2 vectors T_LM|psi> are mutually orthogonal with nonzero squared norm 1/(n+1). They cannot fit in the (n+1)-dimensional symmetric Hilbert space when (q+1)^2>n+1. For n=12,q=3 this would require sixteen orthogonal nonzero vectors in dimension thirteen.

## Exact n=8,t=4 spectral obstruction

Pure symmetric eight-qubit 4-anticoherent states do not exist.

The identity holds for every pure symmetric eight-qubit state. p_J is the nonnegative total-spin-J weight of |psi>⊗|psi>. Four-anticoherence sets a1,...,a4 to zero and would force p0=−4/27. The supplied verifier checks the identity coefficientwise in every allowed total-spin sector using exact fractions.

$$p_0=-\frac4{27}-a_1+\frac{11}{21}a_2+\frac{11}{14}a_3+\frac{143}{126}a_4$$

## Exact n=10,t=5 spectral obstruction

Pure symmetric ten-qubit 5-anticoherent states do not exist.

Five-anticoherence sets a1,...,a5 to zero, while the exact identity would force the nonnegative two-copy weight p6=−27/187. The verifier checks the identity coefficientwise with exact rational arithmetic.

$$p_6=-\frac{27}{187}-\frac{519}{170}a_1-\frac{143}{170}a_2+\frac{754}{255}a_3+\frac{52}{17}a_4+\frac{520}{561}a_5$$

## Exact higher-copy exclusions

The reported exclusions (n,t)=(9,3) and (10,4) are now proved for arbitrary complex pure symmetric states.

For n=9,t=3, twirling four copies yields a positive rotation-invariant operator with unit trace. Eighty-five exact necessary equations leave a three-dimensional linear space of candidate block entries. Four nonnegative diagonal entries obey a strictly positive linear relation; their vanishing forces every parameter to zero. This contradicts unit trace. Two verifiers independently reconstruct the equations and exact ranks. See `ac-exclusions/n9-t3/proof.md`.

For n=10,t=4, an exact degree-six identity expresses the cube of the state norm squared as a sum of twenty cubic multipole invariants, each containing a rank between one and four. Anticoherence makes every summand zero, contradicting a nonzero norm. An independent verifier checks all polynomial coefficients with exact rational arithmetic. See `ac-exclusions/n10-t4/proof.md`.

## Exact n=11,t=4 exclusion and complete core

A hypothetical pure fourth-order eleven-qubit state gives a nonzero positive SU(2)-invariant three-copy moment on Sym^3(C^12). Its 27 real multiplicity coordinates obey 27 necessary rational equations. The supplied integer matrix has determinant 3 modulo 37, hence is invertible. It forces the moment to vanish, a contradiction. Two verifiers independently reconstruct the necessary equations. Taking real parts of positive blocks preserves the equations and nonzero trace, covering arbitrary complex states.

A second certificate is a 26-term degree-six polynomial identity whose left-hand side is the cube of the squared state norm and whose right-hand terms all contain a multipole of rank one through four. An independent rational polynomial expansion verifies all 5,664 monomial coefficients. See `ac-exclusions/n11-t4/proof.md` and `ac-exclusions/n11-t4/moment/proof.md`.

Highest verified orders for n=2,...,12 are 1,1,2,1,3,2,3,2,3,3,5, and all are now proved maximal. The small-system catalogue is complete. In particular, pure fourth-order AC first occurs at n=12.

The earlier two-copy relaxation points remain mathematically feasible. The stronger higher-copy certificates exclude physical states at (9,3), (10,4), and (11,4); feasibility of a relaxation does not imply existence.

Run `python3 verify_small_n_audit.py` for the core witnesses and two-copy identities. All higher-copy proof packages provide a neighboring `verify.py`; the N9 and N11 moment packages also include independent verifiers. Python 3 and SymPy suffice.
