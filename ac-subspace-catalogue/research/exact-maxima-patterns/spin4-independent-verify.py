"""Independent Cartesian-moment audit of the auxiliary-spin-9/2 bound."""
from pathlib import Path
import json
import sympy as sp

j, x, z, q, Q = sp.symbols("j x z q Q")
A, B = sp.symbols("A B")
solution = sp.solve([12*A+3*B-q*q, 6*A+9*B-q*(q-1)], [A, B])
assert sp.factor(solution[A]-q*(2*q+1)/30) == 0
assert sp.factor(solution[B]-q*(q-2)/15) == 0
quartic = sp.factor(2*solution[A]*Q*Q+solution[B]*Q*(Q-1))
assert sp.factor(quartic-q*Q*((3*q-1)*Q-q+2)/15) == 0
generic_moments = [sp.S.One, sp.S.Zero, q*Q/3, -q*Q/6, quartic]
moments = [sp.factor(m.subs({q:j*(j+1), Q:sp.Rational(99,4)}))
           for m in generic_moments]

s = sp.Rational(9,2)
nodes = [j*r+(r*(r+1)-s*(s+1))/2 for r in [-s+i for i in range(10)]]
root_indices = [0,3,4,8]
polynomial = sp.Poly(sp.prod(x-nodes[i] for i in root_indices), x)
trace_from_cartesian_moments = sp.factor(sum(
    polynomial.nth(a)*moments[a] for a in range(5)))
claimed_trace = 9*(j+1)*(4*j**3+724*j*j+2637*j+6750)/20
assert sp.factor(trace_from_cartesian_moments-claimed_trace) == 0
top = sp.factor(polynomial.as_expr().subs(x,nodes[-1])-claimed_trace)
claimed_top = 9*j*(596*j**3+4972*j*j+14489*j+11688)/20
assert sp.factor(top-claimed_top) == 0

for b in [claimed_trace, top]:
    assert all(c>0 for c in sp.Poly(b,j).all_coeffs() if c != 0)
signs = []
for index, node in enumerate(nodes):
    value = sp.factor(polynomial.as_expr().subs(x,node)-claimed_trace)
    positive = value if index == 9 else -value
    coefficients = sp.Poly(sp.expand(positive.subs(j,z+s)),z).all_coeffs()
    assert all(c>0 for c in coefficients)
    signs.append("positive" if index == 9 else "negative")

report = {
    "independent_audit_passed": True,
    "theorem": "10k <= N+10 for AC order >=4 and N>=9",
    "trace_derivation": "Cartesian contractions and Casimir commutators; no total-spin trace sum",
    "fourth_moment_generic": str(quartic),
    "scalar_trace": str(trace_from_cartesian_moments),
    "top_eigenvalue_after_scalar_subtraction": str(top),
    "spectral_signs_bottom_to_top": signs,
    "positive_sector_dimension": "N+10",
    "all_signs_strict_for_entire_domain": True,
    "scope": "Upper bound, with no general attainment claim"
}
payload = json.dumps(report,indent=2)+"\n"
Path(__file__).with_name("spin4-independent-verification.json").write_text(payload)
print(payload,end="")
