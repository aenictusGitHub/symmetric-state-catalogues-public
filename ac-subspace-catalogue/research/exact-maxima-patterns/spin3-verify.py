"""Exact independent checks for the auxiliary-spin-three upper bound."""
from pathlib import Path
import json
import sympy as s
j,x=s.symbols('j x',positive=True)
F=(x+2*j+5)*(x+j+6)*(x-2*j+3)
C=15*(2*j*j+3*j+6)
W=s.Poly(s.expand(F-C),x)
assert s.factor(W.nth(0)+4*j*(j+1)*W.nth(2)-2*j*(j+1)*W.nth(3))==0
expected=[-5*(j+1)*(2*j*j-3*j+18),-C,-C,-2*(2*j+15)*(j*j+j+3),-3*(2*j**3+15*j*j+17*j+30),-C,5*j*(4*j*j+16*j+27)]
for r,value in zip(range(-3,4),expected):
 eigen=j*r+s.Rational(r*(r+1),2)-6
 assert s.factor(W.as_expr().subs(x,eigen)-value)==0
 # j=z+3 covers the theorem domain; all coefficients have the claimed sign.
 z=s.symbols('z',nonnegative=True)
 poly=s.Poly(s.expand((value if r==3 else -value).subs(j,z+3)),z)
 assert all(c>0 for c in poly.all_coeffs())

def spin(J):
 d=int(2*J+1);zz=s.diag(*[J-a for a in range(d)])
 pp=s.zeros(d)
 for a in range(1,d):pp[a-1,a]=s.sqrt(a*(2*J-a+1))
 mm=pp.T
 return ((pp+mm)/2,(pp-mm)/(2*s.I),zz)
S=spin(s.Integer(3));checks=[]
for n in [6,7,8,10]:
 J=s.Rational(n,2);Js=spin(J);d=n+1
 # Physical spin construction, independent of the total-spin eigenvalue formula.
 X=sum((s.kronecker_product(a,b)for a,b in zip(Js,S)),s.zeros(7*d))
 X2=X*X;X3=X2*X
 for power,target in [(X,s.zeros(7)),(X2,4*J*(J+1)*s.eye(7)),(X3,-2*J*(J+1)*s.eye(7))]:
  partial=sum((power[7*a:7*a+7,7*a:7*a+7]for a in range(d)),s.zeros(7))/d
  assert (partial-target).applyfunc(s.simplify)==s.zeros(7)
 checks.append({'n':n,'normalizedPhysicalPartialTracesOfX_X2_X3':True})
result={'passed':True,'theorem':'7k<=N+7 for N>=6 and AC order at least three','symbolicSpectrumChecked':True,'allDomainSignsCertifiedByPositivePolynomialCoefficients':True,'physicalChecks':checks,'scope':'Universal bound, not an existence or optimality statement.'}
Path(__file__).with_name('spin3-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
