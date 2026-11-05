"""Exact checks for the fourth-order auxiliary-spin-9/2 bound."""
from pathlib import Path
import json
import sympy as sp
j,x,z=sp.symbols('j x z',positive=True)
s=sp.Rational(9,2)
rs=[-s+a for a in range(10)]
eigen=lambda r:j*r+(r*(r+1)-s*(s+1))/2
roots=[-s,sp.Rational(-3,2),sp.Rational(-1,2),sp.Rational(7,2)]
F=sp.prod(x-eigen(r) for r in roots)
C=9*(j+1)*(4*j**3+724*j*j+2637*j+6750)/20
mean=sp.factor(sum((2*j+2*r+1)*F.subs(x,eigen(r))for r in rs)/(10*(2*j+1)))
assert sp.factor(mean-C)==0
moments=[1,0,33*j*(j+1)/4,-33*j*(j+1)/8,33*j*(j+1)*(293*j*j+293*j-91)/80]
for a,m in enumerate(moments):
 assert sp.factor(sum((2*j+2*r+1)*eigen(r)**a for r in rs)/(10*(2*j+1))-m)==0
for r in rs:
 w=sp.factor(F.subs(x,eigen(r))-C)
 positive=w if r==s else -w
 # Positive coefficients after j=z+9/2 prove the sign for every j>=9/2.
 assert all(c>0 for c in sp.Poly(sp.expand(positive.subs(j,z+s)),z).all_coeffs())
assert sp.factor(F.subs(x,eigen(s))-C-9*j*(596*j**3+4972*j*j+14489*j+11688)/20)==0

def spin(J):
 d=int(2*J+1);zz=sp.diag(*[J-a for a in range(d)]);pp=sp.zeros(d)
 for a in range(1,d):pp[a-1,a]=sp.sqrt(a*(2*J-a+1))
 mm=pp.T
 return ((pp+mm)/2,(pp-mm)/(2*sp.I),zz)
S=spin(s);checks=[]
for n in [9,10,12]:
 J=sp.Rational(n,2);d=n+1
 X=sum((sp.kronecker_product(a,b)for a,b in zip(spin(J),S)),sp.zeros(10*d))
 power=sp.eye(10*d)
 for a in range(5):
  partial=sum((power[10*b:10*b+10,10*b:10*b+10]for b in range(d)),sp.zeros(10))/d
  target=sp.sympify(moments[a]).subs(j,J)*sp.eye(10)
  assert (partial-target).applyfunc(sp.simplify)==sp.zeros(10)
  if a<4:power=power*X
 checks.append({'n':n,'normalizedPhysicalPartialTracesThroughDegree4':True})
result={'passed':True,'theorem':'10k<=N+10 for N>=9 and AC order at least four','symbolicMeanAndSpectrumChecked':True,'allDomainSignsCertifiedByPositivePolynomialCoefficients':True,'physicalChecks':checks,'scope':'Universal bound, not an existence or optimality statement.'}
Path(__file__).with_name('spin4-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
