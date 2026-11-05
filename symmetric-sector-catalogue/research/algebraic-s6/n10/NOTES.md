# Ten-qubit candidate minimizing the sixth sector

The previously numerical state on Dicke support 0,3,6,9 admits an exact algebraic description. The construction below is derived from the stationary equations, rather than fitted to decimal coefficients. It establishes exact stationarity on the full complex symmetric-state space; it does not establish global minimality.

## Exact normalized state

Write

$$\lvert\psi_{10}\rangle=\frac{\lvert D_{10}^{(0)}\rangle+\sqrt{120u}\lvert D_{10}^{(3)}\rangle+\sqrt{210}\,t\lvert D_{10}^{(6)}\rangle+\sqrt{10u}\,z\lvert D_{10}^{(9)}\rangle}{\sqrt{Q}},\qquad Q=1+120u+210t^2+10uz^2.$$

Here t is the unique real root in the rational interval

`0.34031561845452539230789192544421338373007474986524 < t < 0.34031561845452539230789192544421338373007474986525`

of the irreducible polynomial

$$P_{14}(t)=15693020082432000 t^{14} - 10361344778601600 t^{13} + 2592262933386080 t^{12} - 745743807627504 t^{11} + 279020754020088 t^{10} - 51700201674700 t^{9} + 236498980536 t^{8} + 1359311183406 t^{7} - 181040910612 t^{6} + 4383881496 t^{5} + 687245181 t^{4} - 47112127 t^{3} + 518295 t^{2} + 15399 t + 62.$$

Define

$$z=\frac{A(t)}{B(t)},\qquad u=\frac{30t^2z+8t^2-2tz+20t-1}{8(-30tz+150t+2z^2-29)},$$

where

$$A(t)=389787552000000 t^{14} - 293178516120000 t^{13} + 170727461628000 t^{12} - 15988756474600 t^{11} - 10296253063620 t^{10} + 798022627454 t^{9} + 400513291599 t^{8} - 24502600338 t^{7} - 6495373869 t^{6} + 599281848 t^{5} + 13514649 t^{4} - 2951850 t^{3} + 90026 t^{2} - 249 t - 10,$$

$$B(t)=77957510400000 t^{14} + 13323655926000 t^{13} - 22332190727400 t^{12} + 11442067763140 t^{11} + 236370464214 t^{10} - 555648843266 t^{9} - 15047839098 t^{8} + 8629776678 t^{7} + 341065278 t^{6} - 65621784 t^{5} - 2970030 t^{4} + 382752 t^{3} - 4232 t^{2} - 135 t - 2.$$

Rational interval arithmetic verifies that both recovery denominators are nonzero and u,z are positive throughout the selected root interval. These formulas therefore specify unambiguous real coefficients.

The attained value is exactly

$$S_6=\frac{210(1+168u+66t^2-2uz^2-480ut+96utz)}{Q}.$$

## Decimal values

In increasing excitation order 0,3,6,9, the nonzero normalized amplitudes are

```
0.1503232113481845548317957657843214309324071649491126365454108008198750547831166744219085029237501306
0.5933365153020286287477534248093765184634036898546977493048116331155254865498579877910566139647682141
0.7413402385533004237005004973821598039930247968482861099014826016974413465866704108553409228115281797
0.2752623520224845518811734908572634496650442714767775931412685678720976691107230773504571576337721634
```

They give

`S6 = 73.04625181143282862655926527676680932516076615331620794697006579526721785409907942615617514944364591`.

## Derivation and exact stationarity

Use binomial-scaled variables

`c0=a, c3=√120 b, c6=√210 c, c9=√10 d`.

The norm is `Q=a²+120b²+210c²+10d²`. The sector length is `P/Q²`, where the exact homogeneous quartic P is stored in `polynomials.py` and `system.json`. It follows directly by writing each reduced density matrix in the Dicke basis and applying the purity inversion formula.

Set a=1. The stationary conditions are `a ∂P/∂x - w_x x ∂P/∂a=0` for x=b,c,d and weights 120,210,10. Their primitive integer polynomials are:

```
1200*b**3*c - 232*b**3 - 240*b**2*c*d - 8*b*c**2 - 20*b*c + 16*b*d**2 + b - 30*c**2*d + 2*c*d
2100*b**2*c**2 - 578*b**2*c - 10*b**2 - 420*b*c**2*d - 60*b*c*d + 2*b*d + 168*c**3 + 32*c*d**2 - 3*c
400*b**2*c*d - 96*b**2*d - 120*b*c**2 - 80*b*c*d**2 + 8*b*c + 38*c**2*d + 12*d**3 - d
```

Writing u=b² and z=d/b turns the equations into three polynomials linear in u. Elimination produces P14(c); subresultants recover z=A(c)/B(c), and the first equation recovers u. The verifier checks all three stationary equations modulo P14 using exact rational polynomial arithmetic. It also checks irreducibility, isolates the chosen root by Sturm counting, and certifies the positivity conditions by rational interval arithmetic.

Stationarity in this four-real-coordinate ansatz implies full complex stationarity here. The state is fixed, up to an irrelevant global phase, by collective rotation through 2π/3 about z. Rotation invariance of S6 forces its gradient into the same C3 character, so all missing Dicke components have zero gradient. Complex conjugation leaves S6 invariant, so every imaginary-coordinate derivative vanishes at a real representative. The three ratio equations, together with normalization and phase invariance, therefore cover all remaining variations.

An independent 100-digit matrix calculation gives full normalized 22-real-coordinate gradient norm `1.54591271895748302807552458967e-96`. A numerical full Hessian has five zero symmetry/radial modes and seventeen positive modes, the smallest approximately 17.70836. This supports local minimality, but is not an exact Hessian certificate or a global bound.

## Reproducibility

Copy `polynomials.py`, `algebraic.py`, and `verify.py` into the same directory and run:

```
python3 verify.py
```

SymPy and mpmath are required. Verification regenerates `verified.json` beside the scripts. No external files, pickles, optimizer output, or project-folder paths are needed.

The machine-readable presentation is `parametrization.json`. The scripts `derive.py`, `eliminate.py`, and `parametrize.py` record the derivation and are optional for verifying the final result. The exact value audit below strengthens this: the sector value itself has Galois group S14 and cannot be expressed by radicals over the rationals.


## Exact sector value and obstruction to radical expressions (3 November 2026)

`energy-polynomial.json` gives the monic degree-14 polynomial for the attained
value. Clearing denominators gives a primitive integer polynomial R14(v).
The chosen value is its unique root in

```
73.04625181143282862655926527676680932516076615331620794697
 < S6 <
73.04625181143282862655926527676680932516076615331620794698.
```

`verify_value.py` proves the polynomial identity in the exact field Q[t]/(P14),
checks it against the quartic P/Q², and maps a refined rational interval for t
strictly inside the displayed value interval. It also checks the identity

$$P-\frac{Q P_a}{4}=10080bF_1+10080cF_2+2520dF_3\qquad(a=1).$$

Here P_a denotes the partial derivative with respect to a. This explains why
the shorter rational formula for S6 holds on the stationary branch.

The value polynomial has squarefree reductions with irreducible factor degrees
14 modulo 59, (1,2,11) modulo 83, and (1,13) modulo 151. The first proves
irreducibility over Q. The last two supply a 13-cycle and a transposition in
the Galois group. In degree 14 the 13-cycle excludes nontrivial block systems:
there would have to be 2 or 7 blocks, each of size 7 or 2, all fixed by an
element of order 13, which would then have to fix every point. Thus the group
is primitive. The graph of conjugate transpositions is connected (its
components are blocks), and its edge transpositions generate S14.

The Galois group is consequently S14, which is not solvable. No root of this
irreducible polynomial can be expressed by radicals over Q, even with complex
intermediate radicals. This conclusion concerns the exact attained value;
it does **not** prove that this value is the global minimum.

The standard finite-field cycle criterion and radical-solvability theorem are
reviewed in [Evan Dummit's Galois theory notes, sections 4.4.5–4.4.6](https://dummit.cos.northeastern.edu/docs/fieldthy_4_galois_theory.pdf#page=42).
The complete finite-field factors and interval certificate are in
`value-verification.json`. Reproduce them with:

```
python3 verify_value.py
```
