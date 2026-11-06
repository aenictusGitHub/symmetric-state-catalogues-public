"""Exact symbolic checks for the retained paired/rank-four crease exclusion."""
from pathlib import Path
import sympy as S
import json
v=S.symbols('v',real=True)
q=S.Rational(169,448)*(1-25*v*v)
delta=40*v*v
g=S.Rational(2,7)*(1-5*v)*(135*v-1)
a=S.Rational(1,4)+v/2;h=S.Rational(1,4)+v+5*v*v
k=1+(2*S.sqrt(2)-8)*a+(22-4*S.sqrt(2))*a*a
x=S.Rational(3,2)-S.Rational(474,169)*q
y=S.Rational(99,10)+S.Rational(3234,845)*q
z=11-S.Rational(2300,169)*q
checks={
 'sharp_plane_identity':S.simplify(65*y+63*z-217*x-1011)==0,
 'delta_identity':S.simplify(15+5*x-y-z-delta)==0,
 'g_identity':S.simplify(37+3*x-5*y+48*q+40*v-g)==0,
 'q_omega_endpoint':S.simplify(q.subs(v,S.Rational(3,55))-S.Rational(169,484))==0,
 'q_C5_endpoint':S.simplify(q.subs(v,S.Rational(1,5)))==0,
 'gap_factorization':S.simplify(g*h-4*delta*q-(5*v-1)**2*(305*v*v+126*v-1)/14)==0,
 'k_h_factorization':S.simplify(k-h-(2*v+1)*(1+2*S.sqrt(2)+(2-4*S.sqrt(2))*v)/8)==0,
 'positive_polynomial_at_lower_endpoint':S.simplify((305*v*v+126*v-1).subs(v,S.Rational(3,55))-S.Rational(4102,605))==0,
 'positive_bracket_at_upper_endpoint':S.simplify((1+2*S.sqrt(2)+(2-4*S.sqrt(2))*v).subs(v,S.Rational(1,5))-(7+6*S.sqrt(2))/5)==0}
assert all(checks.values()),checks
report={'status':'PASS','method':'Exact SymPy rational/radical identities; interval inequalities are proved in accompanying text',
        'checks':checks,'proofNote':'audit/N6_boundary_sample_proofs.txt',
        'conclusion':'The actual retained paired-spectral/rank-four curved crease has no attaining state except the C5 junction.'}
p=Path(__file__).with_name('n6-retained-crease-proof-audit.json');p.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
