import sympy as s,json
from sympy.physics.wigner import wigner_6j
from pathlib import Path
x,Q=s.symbols('x Q');c=Q/4;f=c-x/2
lhs=1+3*f/c+30*(f*f-f/2-c*c/3)/(c*(4*c-3))
rhs=3*(10*x*x+(16-12*Q)*x+3*Q*(Q-3))/(Q*(Q-3))
assert s.factor(lhs-rhs)==0
checks=0
for n in range(4,25):
 j=s.Rational(n,2);q=n*(n+2)
 for L in range(n+1):
  exact=sum((-1)**(n+ell+L)*(2*ell+1)*wigner_6j(j,j,ell,j,j,L)for ell in range(3))
  assert exact==rhs.subs({Q:q,x:L*(L+1)})/(n+1);checks+=1
cases=0
for t in range(1,31):
 for k in range(1,21):
  demands=[]
  for h in range(t+1):
   # Twice the tensor rank runs h,h+2,...,t.
   D=sum(ell2+1 for ell2 in range(h,t+1,2));demands.append(k*D+h-1)
  r=t//2
  expected=k*(r+1)*(r+2)if t%2 else (r+1)**2 if k==1 else k*(r+1)**2-1
  assert max(demands)==expected;cases+=1
report={'passed':True,'symbolicFourthOrderPolynomialIdentity':True,'exactSixJChecks':checks,'finiteEnvelopeChecks':cases,'scope':'The accompanying proof establishes all N,t,k; finite checks only cross-check the algebra.'}
Path(__file__).with_name('general-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
