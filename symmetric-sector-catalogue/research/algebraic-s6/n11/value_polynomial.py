import sympy as s,json
from pathlib import Path
R=Path(__file__).resolve().parent;A,B,T,V=s.symbols('A B T V',real=True)
d=json.loads((R/'exact-algebraic.json').read_text());loc={'A':A,'B':B,'T':T}
p=s.sympify(d['Bpolynomial'],locals=loc)
def rem(x):return s.rem(x,p,B)
def inv(x):return s.invert(x,p,B)
ar=rem(s.sympify(d['AnumeratorReduced'],locals=loc)*inv(s.sympify(d['AdenominatorReduced'],locals=loc)))
def evalpoly(poly):
 poly=s.Poly(poly,A,B,T);rr=0
 for (ia,ib,it),coef in poly.terms():rr=rem(rr+coef*rem(ar**ia)*B**ib*rem(tr**it))
 return rr
tn=s.sympify(d['Tnumerator'],locals=loc);td=s.sympify(d['Tdenominator'],locals=loc)
tr=rem(rem(tn.subs(A,ar))*inv(rem(td.subs(A,ar))))
Q=evalpoly(A*A*T+165*B*B+462*T+55)
F=462*A**4*T**2+152460*A*A*T*B*B+49896*A*A*T*T-4620*A*A*T-1774080*A*B*B*T+1774080*A*B*T+7477470*B**4+16659720*B*B*T+9438660*B*B-39029760*B*T+77810040*T*T+8842680*T+732270
vr=rem(evalpoly(F)*inv(rem(Q*Q)))
cols=[]
for i in range(11):
 poly=s.Poly(rem(vr*B**i),B);cols.append(s.Matrix([poly.nth(j) for j in range(11)]))
M=s.Matrix.hstack(*cols)
ch=M.charpoly('V')
char=s.Poly(ch.as_expr(),ch.gen).clear_denoms()[1].primitive()[1]
print(char.as_expr(),flush=True)
(R/'value-polynomial.json').write_text(json.dumps({'variable':'V','polynomial':str(char.as_expr()),'irreducible':bool(char.is_irreducible)},indent=2)+'\n')
