#!/usr/bin/env python3
"""Exact, portable audit of the N=14, k=3, t=2 radical construction.

Every ordered cross-reduction entry is checked directly by Dicke splitting.
The verifier uses the quadratic field for c and the two real square roots f,g;
all reductions are rational polynomial arithmetic, without floating-point tests.
Requires SymPy. Reads only basis.json next to this file.
"""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent

def verify():
 raw=(ROOT/'basis.json').read_bytes();record=json.loads(raw)['subspaces'][0]
 c,f,g=s.symbols('c f g',real=True);env={'c':c,'f':f,'g':g}
 P=s.sympify(record['parameterRelations']['cPolynomial'],locals=env)
 u=s.sympify(record['parameterRelations']['fSquared'],locals=env)
 v=s.sympify(record['parameterRelations']['gSquared'],locals=env)
 norm=s.sympify(record['parameters']['normalizationSquared'],locals=env)
 c0=s.sympify(record['parameters']['c']);f0=s.sqrt(u.subs(c,c0));g0=s.sqrt(v.subs(c,c0));norm0=norm.subs(c,c0)
 assert s.simplify(P.subs(c,c0))==0
 lo,hi=map(s.Rational,record['parameterRelations']['cIsolatingInterval'])
 assert s.Poly(P,c).count_roots(lo,hi)==1 and lo<c0<hi
 # These are affine functions, so positivity throughout the isolating interval
 # follows from exact positivity at its two endpoints.
 for expr in [u,v,norm]:
  assert s.Poly(expr,c).degree()<=1
  assert expr.subs(c,lo)>0 and expr.subs(c,hi)>0
 def field_reduce(expr):
  num,den=s.fraction(s.cancel(expr));poly=s.Poly(num,f,g)
  reduced=s.S.Zero
  for (i,j),coef in poly.terms():reduced+=coef*u**(i//2)*v**(j//2)*f**(i%2)*g**(j%2)
  parts={}
  for powers,coef in s.Poly(s.expand(reduced),f,g).terms():
   numerator,denominator=s.fraction(s.cancel(coef))
   remainder=s.rem(numerator*s.invert(denominator,P,c),P,c)
   if remainder!=0:parts[str(powers)]=str(s.factor(remainder))
  # Every denominator used by the input is a power of f times a nonzero
  # rational constant. c-dependent denominators introduced by reduction are
  # inverted only after exact gcd verification via invert().
  if den.has(f) or den.has(g):
   mon=s.Poly(den,f,g);assert len(mon.terms())==1
   for val in mon.coeffs():assert val!=0
  return parts
 rows=[];coefficients_checked=0
 for basis in record['basis']:
  row=[s.sympify(a,locals=env)for a in basis['weightedCoefficientExpressions']]
  assert len(row)==15
  for r,(weighted,saved)in enumerate(zip(row,basis['exactCoefficients'])):
   expected=s.sqrt(s.binomial(14,r))*weighted.subs({c:c0,f:f0,g:g0})/s.sqrt(norm0)
   assert s.simplify(s.sympify(saved)-expected)==0,(basis['id'],r,'closed-radical mismatch')
   coefficients_checked+=1
  rows.append(row)
 gram_count=0
 for a,A in enumerate(rows):
  for b,B in enumerate(rows):
   expr=sum(s.binomial(14,r)*A[r]*s.conjugate(B[r])for r in range(15))-int(a==b)*norm
   assert not field_reduce(expr),(a,b,'Gram');gram_count+=1
 count=0
 for a,A in enumerate(rows):
  for b,B in enumerate(rows):
   for r in range(3):
    for q in range(3):
     # The omitted factor sqrt(C(2,r) C(2,q))/norm is nonzero.
     expr=sum(s.binomial(12,l)*A[r+l]*s.conjugate(B[q+l])for l in range(13))
     if a==b and r==q:expr-=norm/(3*s.binomial(2,r))
     assert not field_reduce(expr),(a,b,r,q,'cross reduction')
     count+=1
 next_failure=None
 for a,A in enumerate(rows):
  for b,B in enumerate(rows):
   for r in range(4):
    for q in range(4):
     expr=sum(s.binomial(11,l)*A[r+l]*s.conjugate(B[q+l])for l in range(12))
     if a==b and r==q:expr-=norm/(4*s.binomial(3,r))
     residual=field_reduce(expr)
     if residual and next_failure is None:next_failure={'basisA':a,'basisB':b,'row':r,'column':q,'nonzeroReducedFieldNumerator':residual}
 assert next_failure is not None
 return {'passed':True,'id':record['id'],'basisSha256':hashlib.sha256(raw).hexdigest(),'method':'Every ordered Dicke cross reduction checked by exact rational polynomial reduction modulo the quadratic for c and f²=u(c), g²=v(c); saved closed radicals checked against that parameterization.','closedRadicalCoefficientsChecked':coefficients_checked,'orthonormalGramEntriesChecked':gram_count,'fullCrossReductionEntriesChecked':count,'achievedOrder':2,'nextOrderFails':True,'nextOrderCounterexample':next_failure,'positiveRadicandsCertifiedOnIsolatingInterval':True,'dimension':3,'maximumDimension':3,'upperBound':'Tensor-Gram: k*(floor(t/2)+1)^2 <= N+1 gives k<=floor(15/4)=3.'}

if __name__=='__main__':
 report=verify();(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
