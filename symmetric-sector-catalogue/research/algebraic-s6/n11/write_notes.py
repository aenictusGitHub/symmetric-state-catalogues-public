import json,sympy as s
from pathlib import Path
R=Path(__file__).resolve().parent;d=json.loads((R/'exact-algebraic.json').read_text());v=json.loads((R/'verification.json').read_text())
B=s.symbols('B',real=True)
def tex(key):return s.latex(s.sympify(d[key],locals={'B':B}))
text=r'''Exact algebraic Dicke coefficients for the N=11 minimum-S6 candidate
==================================================================

This construction defines the numerically best-known candidate exactly.
It proves attainment of its algebraic sector value and stationarity in the
full complex symmetric-state space. It does not prove global optimality.

Let B be the unique real root in the rational interval

  0.2498461150300332175000977146562933123120
    < B <
  0.2498461150300332175000977146562933123121

of the following irreducible degree-11 polynomial:

P(B) = '''+d['Bpolynomial']+r'''.

Define A = N(B)/D(B), with

N(B) = '''+d['AnumeratorReduced']+r''',

D(B) = '''+d['AdenominatorReduced']+r'''.

Now set

T = [10 A B^2 - 10 A B + 706 B^2 + 440 B + 39]
    / [A^2 - 84 A B + 1848 B + 917],

Q = A^2 T + 165 B^2 + 462 T + 55.

The state, in the normalized Dicke basis, is exactly

|psi_11> = [ A sqrt(T) |D_11^(0)> + sqrt(165) B |D_11^(3)>
             + sqrt(462 T) |D_11^(6)> + sqrt(55) |D_11^(9)> ] / sqrt(Q).

Here A < 0 and T > 0. All other Dicke coefficients vanish.
The root construction is derived from exact polynomial elimination;
the coefficients were not obtained by fitting decimal values.

100-significant-digit coefficients:
'''
for r,c in v['coefficients100digits'].items():text+=f'  c_{r} = {c}\n'
text+='\nSector value:\n  S_6 = '+v['sectorValue100digits']+'\n'
value_data=json.loads((R/'value-polynomial.json').read_text())
text+='\nEquivalently, this attained value is the unique real root V in\n'
text+='  113.1090435663466909908720436191490453102178 < V\n'
text+='    < 113.1090435663466909908720436191490453102179\n'
text+='of the irreducible degree-11 polynomial\n\n  '+value_data['polynomial']+' = 0.\n'
text+=r'''
Derivation and exact verification
--------------------------------

For real coefficients supported on r=0,3,6,9, rescale
c_r = sqrt(C(11,r)) x_r and fix x_9=1. Write x_0=A sqrt(T),
x_3=B, x_6=sqrt(T). The unnormalized norm is Q. The three
independent constrained-stationarity equations become

  T A^3 - 16 T A^2 B - 19 A B^2 + 352 T A B - 110 T A
    - 31 A - 880 B^2 + 880 B = 0,

  T(3 A^2 B - 48 A B^2 - 32 A B + 16 A + 1056 B^2
    - 178 B - 352) - 241 B^3 + 91 B = 0,

  T(A^2 - 84 A B + 1848 B + 917)
    - 10 A B^2 + 10 A B - 706 B^2 - 440 B - 39 = 0.

Eliminating A,T produces P(B) as the factor selecting this solution.
The rational function for A comes from the linear subresultant and is
reduced modulo P(B). Each of the three stationary equations has exactly
zero remainder in Q[B]/(P). SymPy verifies that P is irreducible and that
the stated interval contains exactly one real root.

The algebraic degree is substantive. Modulo 71, P has two distinct
irreducible factors of degrees 2 and 9:

  (B^2+26B+28)
  (B^9-24B^8-28B^7+5B^6+33B^5-3B^4+15B^3+B^2-6B+35),

up to a nonzero constant factor. Hence its Galois group contains a
permutation with cycle lengths 2 and 9; its ninth power is a transposition.
Irreducibility gives a transitive group of prime degree 11, which is
primitive. A primitive permutation group containing a transposition is
the full symmetric group, so Gal(P)=S_11. Thus B is not expressible by
radicals over Q. In particular, these four nonzero Dicke coefficients
cannot all have radical expressions, because B=c_3/(sqrt(3)c_9).
The exact isolated-root representation is therefore appropriate.

The full complex stationarity follows as well. The state is invariant,
up to a global phase, under the collective rotation through 2 pi/3 about z.
The gradient can therefore have support only on Dicke indices 0 mod 3.
For a real state, imaginary-direction derivatives vanish by complex
conjugation invariance. The remaining real derivatives are precisely the
three equations above after normalization. An independent 110-digit
evaluation of the full projected gradient gives maximum residual
'''+v['fullRealGradientResidual']+r'''.

Exact family from time reversal
------------------------------

Define the antiunitary time-reversal map by

  Theta |D_n^(r)> = (-1)^r |D_n^(n-r)>,

with complex conjugation of scalar coefficients. For odd n,
Theta^2=-I, so |psi> and Theta|psi> are orthonormal. Every even-weight
Pauli product P obeys Theta P Theta^{-1}=P. Its diagonal matrix elements
in that two-dimensional span coincide, and its off-diagonal elements
vanish. For example the skew relation

  <x|Theta y> = -<y|Theta x>

with y=P x proves <psi|P Theta psi>=0. Therefore any normalized state

  |psi(alpha,beta)> = alpha |psi> + beta Theta|psi>,
  |alpha|^2 + |beta|^2 = 1,

has the same individual even-weight Pauli expectation values and hence
the same S_2,S_4,S_6,S_8,S_10. This gives two additional real flat
directions after removal of the overall phase. For the state above,

Theta|psi_11> = [ -sqrt(55) |D_11^(2)> + sqrt(462 T) |D_11^(5)>
                 -sqrt(165) B |D_11^(8)> + A sqrt(T) |D_11^(11)> ] / sqrt(Q).

Files and reproducibility
-------------------------

derive.py                 Exact quartic sector polynomial and high precision root.
eliminate.py              Exact projective stationary ideal and Groebner basis.
reduce.py                 Reduced equations and resultant.
rational_univariate.py    Degree-11 root and rational A(B).
verify_exact.py            Root isolation, exact quotient-field checks, full gradient.
exact-algebraic.json       Machine-readable exact expressions.
verification.json          High precision coefficients and verification results.
value-polynomial.json      Optional degree-11 polynomial for the attained S6 value.

Run with Python 3, SymPy and mpmath. The optional family check also needs NumPy.
'''
(R/'analytic-notes.txt').write_text(text)
latex=r'''\begin{align}
|\psi_{11}\rangle&=\frac{A\sqrt{T}|D_{11}^{(0)}\rangle+\sqrt{165}B|D_{11}^{(3)}\rangle+\sqrt{462T}|D_{11}^{(6)}\rangle+\sqrt{55}|D_{11}^{(9)}\rangle}{\sqrt{Q}},\\
T&=\frac{10AB^2-10AB+706B^2+440B+39}{A^2-84AB+1848B+917},\\
Q&=A^2T+165B^2+462T+55,\qquad A=\frac{N(B)}{D(B)}.
\end{align}
'''
latex+='\nP(B)='+tex('Bpolynomial')+'=0.\n\nN(B)='+tex('AnumeratorReduced')+'.\n\nD(B)='+tex('AdenominatorReduced')+'.\n'
(R/'formula.tex').write_text(latex)
print('Wrote analytic-notes.txt and formula.tex')
