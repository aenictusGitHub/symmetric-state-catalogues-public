# A rigorous obstruction when the four doubled points form antipodal pairs

The result proved here is restricted: a pure symmetric 12-qubit state of full
Majorana partition `(2,2,2,2,1,1,1,1)` cannot be 3-anticoherent if its four doubled
points form two antipodal pairs. The four simple points are unrestricted. This
**does not** settle the partition for arbitrary doubled-point geometry.

The proof uses exact rational polynomial identities and real-root counts. Run
`python verify-antipodal-obstruction.py` from this directory (SymPy is required).
The output is `antipodal-obstruction-verification.json`. The verifier reconstructs
its matrices and equations rather than importing the discovery scripts.

## 1. Complete normal form for the stated family

Two distinct unoriented antipodal axes lie in a plane. A spatial rotation sends
that plane to the equator. A further rotation about the polar axis places the
four points at `±exp(iθ), ±exp(−iθ)` in the root coordinate. Distinctness gives
`0<θ<π/2`. Consequently the double-root factor can be written

\[
 Q(z)=z^4+2hz^2+1,\qquad h=-\cos(2\theta)\in(-1,1).
\]

Every state in the stated family therefore has Majorana polynomial
`P=Q²R`, with an arbitrary complex quartic `R`. A degree below four is allowed;
the missing degree records roots at infinity. Common roots of `Q` and `R`,
repeated roots of `R`, and the endpoints `h=±1` cannot occur in the target full
partition. We permit these degenerations during the algebraic calculation and
exclude them when identifying the surviving partition.

We use the normalized-Dicke convention

\[
 P(z)=\sum_{r=0}^{12}(-1)^r\sqrt{\binom{12}{r}}c_rz^r
\]

up to a nonzero scalar. Rotations and global phase do not affect anticoherence
or multiplicities.

## 2. Third-order anticoherence forces the residual quartic to be antipodal

Consider the time-reversal-real quartic basis

\[
 B=(z^4+1,\ z^2,\ i(z^4-1),\ z^3-z,\ i(z^3+z)).
\]

It is a complex basis of all quartics, including polynomials of lower degree.
Write `R=Σ(x_j+i y_j)B_j`, with real vectors `x,y`. Since `Q` is time-reversal
real, the corresponding spin-six vectors `W_j=Q²B_j` are also time-reversal real.
The odd part of the pure density operator is therefore a linear combination of

\[
 i(W_iW_j^\dagger-W_jW_i^\dagger),\quad i<j,
\]

with coefficients given, up to a common sign, by `x_i y_j−y_i x_j`.

Taking a three-qubit reduced state commutes with time reversal. The identity
`ρ₃=I₄/4` is time-reversal even, so the reduced image of this odd part must vanish.
The verifier constructs all 16 real components of this reduction as a `16×10`
matrix `M(h)`. Explicitly, before positive row scalings, the matrix entries follow
from

\[
 \rho_{ac}=(-1)^{a+c}\sqrt{\binom3a\binom3c}
 \sum_{q=0}^9\frac{\binom9q p_{a+q}\overline{p_{c+q}}}
 {\binom{12}{a+q}\binom{12}{c+q}}.
\]

The selected ten rows listed in `odd-restriction.json` have determinant

\[
 -697827328000000(h^2+3)^2 A_0(h)A_-(h)A_+(h),
\]

where the three degree-twelve factors are recorded in that file and in the
verification output. The verifier recomputes the determinant and uses exact
Sturm counts to establish that none of its factors has a real root in `[-1,1]`.
Thus `M` has rank ten throughout the allowed interval. All wedge coordinates
`x_i y_j−y_i x_j` vanish, so `x,y` are proportional. A global phase makes `R`
time-reversal real.

This step does not assume that the four simple points were antipodal. It proves
that they would have to be antipodal in any candidate 3-AC state in this family.

## 3. Exact elimination in the time-reversal-real family

After the phase choice, put

\[
 R=(b+id)z^4+(e+if)z^3+cz^2+(-e+if)z+b-id,
 \qquad b,c,d,e,f\in\mathbb R.
\]

Odd moments vanish automatically. Second-order isotropy is equivalent to the
five real equations from `⟨J_z²⟩=14`, `⟨{J_z,J_+}⟩=0`, and `⟨J_+²⟩=0`.
Their primitive integer polynomials are in `antipodal-reduction.json`; the
verifier independently reconstructs them from the Dicke coefficients.

* If `ef≠0`, set `t=f/e` and solve the two linear off-diagonal equations for
  `b,c` in terms of `d`. The remaining off-diagonal equation is

  \[
  t^2D(h)K(h)e^2+16(h^2+1)
       [t^2A(h)+A(-h)]d^2=0,
  \]

  where
  \[
  D=(h^2-2h+5)(h^2+2h+5),\quad K=7h^4+18h^2+7,
  \]
  \[
  A=h^6+2h^5+9h^4+12h^3+51h^2+34h+115.
  \]
  Here `D,K>0` and `A(±h)≥67` on `[-1,1]`. Since `e≠0`, this is impossible.
  The derivation does not divide by `d`, so it includes `d=0`.

* If exactly one of `e,f` is nonzero, another off-diagonal equation forces `d=0`.
  For example, when `e=0`, its coefficient is
  `h⁴+4h³+12h²+16h+23≥3`. The remaining linear equation fixes
  `c=2b(2h²−h+9)/(h²−h+2)`. The two diagonal conditions then form a homogeneous
  linear system in `b²,f²`. Its determinant is

  \[
  48(h-3)(h+1)^3(h^2+3)(h^2-h+2)^2
       (h^3-4h^2+h-2)(h^3+h^2-h-17).
  \]

  It has no root in `(-1,1)`, verified exactly. The `f=0` branch is also checked
  directly. Both branches contradict the nonzero active component.

* If `e=f=0` but `d≠0`, the last equation fixes
  `c=−2bh/(h²+2)`. The diagonal equation, after multiplying by its positive
  denominator, becomes

  \[
  (-16h^8+1104h^4+4320h^2+4320)b^2
  +(-16h^8-48h^6+1088h^4+4416h^2+4352)d^2=0.
  \]

  The coefficients are at least `4304` and `4288` on `[-1,1]`, respectively.
  This is impossible.

* In the remaining branch `d=e=f=0`, suppose first that `b≠0` and rescale it
  to one. Eliminating `c` gives

  \[
  3072h^2(h-3)^2(h+3)^2(h^2+3)^2
       (8h^6+51h^4+258h^2+195)=0.
  \]

  The only root in `(-1,1)` is `h=0`, where the first equation is `1080=0`.
  Hence `b=0`. The nonzero-state condition now gives `c≠0`; the diagonal
  equation is `−2c²h²(3h²+11)=0`, forcing `h=0`.

## 4. Multiplicity guard and conclusion

The sole surviving projective polynomial is

\[
 P(z)=z^2(z^4+1)^2.
\]

It has four simple roots of `z⁴+1`, each doubled, a double root at zero, and a
double root at infinity because its affine degree is ten while `n=12`. Its full
partition is therefore `(2,2,2,2,2,2)`, the double octahedron. The verifier checks
both this square-free decomposition and `ρ₃=I₄/4` exactly.

The survivor lies outside `(2⁴,1⁴)`. Therefore any hypothetical 3-AC state in the
target partition must have a doubled-point configuration that is **not** two
antipodal pairs. No global nonexistence claim follows from this restriction.
