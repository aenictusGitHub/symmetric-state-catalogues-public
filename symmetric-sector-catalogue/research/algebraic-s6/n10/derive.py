import sympy as s
from math import comb
from pathlib import Path
import json,time
n=10;k=6
a,b,c,d=s.symbols('a b c d');xs={0:a,3:b,6:c,9:d};vs=[a,b,c,d];weights=[comb(n,i) for i in xs]
Q=sum(w*x*x for w,x in zip(weights,vs))
def purity(r):
 total=0
 for i in range(r+1):
  for j in range(r+1):
   elem=sum(comb(n-r,q)*xs.get(i+q,0)*xs.get(j+q,0) for q in range(n-r+1))
   if elem!=0:total+=comb(r,i)*comb(r,j)*elem**2
 return s.expand(total)
P=s.expand(comb(n,k)*sum((-1)**(k-r)*2**r*comb(k,r)*purity(min(r,n-r)) for r in range(k+1)))
print('Q=',Q,flush=True);print('P=',P,flush=True)
F=[s.factor(s.diff(P,x)*a-w*x*s.diff(P,a)).subs(a,1) for x,w in zip(vs[1:],weights[1:])]
F=[s.Poly(f,b,c,d).primitive()[1].as_expr() for f in F]
print('F=',F,flush=True)
guess=[.5933365155/(comb(10,3)**.5)/.1503232039,.7413402194/(comb(10,6)**.5)/.1503232039,.2752624072/(comb(10,9)**.5)/.1503232039]
sol=s.nsolve(F,[b,c,d],guess,prec=100,tol=s.Float('1e-95'),maxsteps=100)
val=s.N(P.subs(a,1).subs(dict(zip([b,c,d],sol)))/Q.subs(a,1).subs(dict(zip([b,c,d],sol)))**2,100)
coeff=[s.N(s.sqrt(w)*x/s.sqrt(Q.subs(a,1).subs(dict(zip([b,c,d],sol)))),100) for w,x in zip(weights,[1,*sol])]
print('ratios=',sol,flush=True);print('value=',val,flush=True);print('coeff=',coeff,flush=True)
result={'n':n,'k':k,'support':[0,3,6,9],'Q':str(Q),'P':str(P),'stationaryPolynomials':[str(f) for f in F],'ratioValues':[str(x) for x in sol],'value':str(val),'coefficients':[str(x) for x in coeff]}
(Path(__file__).resolve().parent / 'system.json').write_text(json.dumps(result,indent=2))
# Save objects in Python-readable source for elimination.
(Path(__file__).resolve().parent / 'polynomials.py').write_text('import sympy as s\na,b,c,d=s.symbols("a b c d")\nQ='+repr(Q)+'\nP='+repr(P)+'\nF='+repr(F)+'\n')
