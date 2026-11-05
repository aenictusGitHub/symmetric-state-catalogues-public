# Nine-qubit candidate maximizing the fifth sector

The numerical four-component candidate has an exact nested-radical representative. Its attained value is

$$S_5=\frac{130928+6832\sqrt{61}}{1461}=126.1380053290857377280590456746791107854\ldots.$$

The coefficient formulas and full complex stationarity are verified exactly. Global maximality remains conjectural.

## Explicit Dicke coefficients

Define

$$t=\frac{\sqrt{2447347+310013\sqrt{61}}-991-159\sqrt{61}}{2050},\qquad z=5(15-\sqrt{61})t,$$

$$u=\frac{t(4tz-29t-z+8)}{-84tz+336t+z^2-25},\qquad Q=1+84u+84t^2+uz^2.$$

Then a normalized representative is

$$|\psi_9\rangle=\frac{|D_9^{(0)}\rangle+\sqrt{84u}|D_9^{(3)}\rangle+\sqrt{84}\,t|D_9^{(6)}\rangle+z\sqrt u|D_9^{(9)}\rangle}{\sqrt Q}.$$

Here `t<0`, `z<0`, and `u,Q>0`. The positive square roots are intended. In increasing excitation order 0,3,6,9, the normalized amplitudes are

```
0.8646511015018027881285720191021845476414517627367070838529107795636248890409833499751085642429418346
0.4913278799569659477595176014260719122101759821215476537421365357643562748453873715808879320901478241
-0.1017966652597417262313889826793553797792328982004191568822266189236706318298007660920772927623067731
-0.02475532247231862510547091825872476732328651365133739283646591094367446366723831468504614059460365234
```

All other Dicke coefficients vanish.

## Exact algebraic verification

Using the binomial-scaled variables `c0=a`, `c3=√84 b`, `c6=√84 c`, `c9=d`, write the sector as `P/(a²+84b²+84c²+d²)²`. The quartic P and its stationary equations are in `polynomials.py`; they are derived directly from the reduced-purity inversion formula in `derive.py`.

Setting `a=1`, `u=b²`, and `z=d/b` makes the three stationary equations linear in u. Eliminating u,z yields a factor

$$P_4(t)=20500t^4+39640t^3-10178t^2+252t+5.$$

This polynomial is irreducible over the rationals. The relevant quadratic factor over the field containing √61 is

$$t^2+\frac{991+159\sqrt{61}}{1025}t+\frac{15+\sqrt{61}}{820}=0.$$

Its larger root is exactly the displayed nested radical. A Sturm count isolates this root in

`−0.01284555902806194669499550925759084638053248768653 < t < −0.01284555902806194669499550925759084638053248768652`.

The verifier checks that the radical lies in this interval using rational square-root interval arithmetic. It verifies both recovery denominators are nonzero and that u,Q are positive. Exact polynomial remainders establish all stationary identities. The attained value obeys

$$1461S_5^2-261856S_5+9784320=0,$$

and the selected branch is the larger root, `(130928+6832√61)/1461`. Equality of the objective with this radical is checked algebraically, not inferred from matching decimals.

## Full-space stationarity and local check

Collective rotation about z through 2π/3 fixes this state up to global phase. Because S5 is rotation invariant, its gradient has the same C3 character, forcing the missing Dicke components to have zero derivative. Complex conjugation leaves S5 invariant; thus all imaginary-coordinate derivatives vanish for this real representative. The exact support equations therefore establish stationarity throughout the full complex symmetric-state space, not only within a numerical four-amplitude restriction.

An independent 100-digit reduced-density-matrix calculation gives normalized full complex gradient norm `7.23701491508115784635308218615e-99`. The numerical Hessian has fifteen negative physical eigenvalues and five zero normalization/phase/rotation modes. The least negative physical eigenvalue is approximately −1.38298935449. These checks support a strict local maximum modulo the symmetries; they are not a proof of global maximality.

Collective SU(2) rotations and overall phase preserve the value. The Kramers superposition construction that preserves even sector lengths does not preserve S5, so it does not generate an equal-value family for this result.

## Reproducibility

Keep `polynomials.py`, `algebraic.py`, `radical.py`, and `verify.py` together, then run

```
python3 verify.py
```

Only SymPy and mpmath are required. All paths are relative to the script directory. `verified.json` records the checks; `parametrization.json` contains the compact formulas and full-precision values for the website. The derivation scripts are optional for verifying the final formulas.
