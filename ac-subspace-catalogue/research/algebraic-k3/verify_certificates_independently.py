#!/usr/bin/env python3
"""Independent exact audit of the N15 and N17 algebraic-root certificates.
No case-specific definition.py or verify.py is imported. Polynomial evaluation,
Jacobians, Hessians, Newton bounds and row reduction use fractions.Fraction.
SymPy is used only to parse saved symbolic coefficient expressions.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parent

def rational(v):return Q(str(v))
def poly(expression,variables):
 return {powers:rational(coef) for powers,coef in s.Poly(s.expand(expression),*variables).terms() if coef}
def evaluate(p,point):
 return sum(coef*product(point[j]**power for j,power in enumerate(powers)if power)for powers,coef in p.items())
def product(values):
 z=Q(1)
 for x in values:z*=x
 return z
def derivative(p,j):
 out={}
 for powers,coef in p.items():
  if powers[j]:
   a=list(powers);a[j]-=1;out[tuple(a)]=coef*powers[j]
 return out
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def maxrows(M):return max(sum(abs(x)for x in row)for row in M)
def echelon(polys):
 pivots=[]
 for source in polys:
  row=source.copy()
  for pivot,basis in pivots:
   if pivot in row:
    factor=row[pivot]
    for key,value in basis.items():
     row[key]=row.get(key,Q(0))-factor*value
     if not row[key]:del row[key]
  if row:
   pivot=min(row);factor=row[pivot];row={k:v/factor for k,v in row.items()};pivots.append((pivot,row))
 return pivots
def in_span(source,pivots):
 row=source.copy()
 for pivot,basis in pivots:
  if pivot in row:
   factor=row[pivot]
   for key,value in basis.items():
    row[key]=row.get(key,Q(0))-factor*value
    if not row[key]:del row[key]
 return not row
def decimal_upper(value):
 if value==0:return '0'
 exponent=0
 while value<Q(1,10**(exponent+1)):exponent+=1
 return '1e-'+str(exponent) if exponent else str(value.numerator//value.denominator+1)

def audit(n):
 folder=ROOT/('n'+str(n));basis_raw=(folder/'basis.json').read_bytes();cert_raw=(folder/'certificate.json').read_bytes()
 basis=json.loads(basis_raw);cert=json.loads(cert_raw);definition=basis['algebraicDefinition']
 variables=s.symbols(' '.join(definition['variables']),real=True);env={str(x):x for x in variables};D=len(variables)
 if n==15:
  assert cert['variables']==list(env)
  equations=[{tuple(t['powers']):rational(t['coefficient'])for t in p}for p in cert['equations']]
  assert equations==[poly(s.sympify(e,locals=env),variables)for e in definition['equations']]
  inverse=cert['inverseApproximation'];fixed={}
 else:
  equations=[poly(s.sympify(e,locals=env),variables)for e in cert['equations']]
  inverse=cert['inverseJacobianApproximation'];fixed=definition['fixedCoordinates']
 assert len(equations)==D and all(sum(powers)<=2 for p in equations for powers in p)
 center=list(map(rational,cert['center']));radius=rational(cert['radius']);R=[list(map(rational,row))for row in inverse]
 assert len(center)==len(R)==D and all(len(row)==D for row in R)
 for name,value in zip(env,center):assert list(map(rational,definition['parameterIntervals'][name]))==[value-radius,value+radius]
 for index,value in fixed.items():assert list(map(rational,definition['parameterIntervals']['x'+index]))==[rational(value)]*2
 # Independently rebuild all full cross conditions from the SAVED exact Dicke
 # coefficients, instead of the selected equations or their generating ansatz.
 vectors=[];render_error=0.0;subs={x:s.Rational(value.numerator,value.denominator)for x,value in zip(variables,center)}
 for record in basis['basis']:
  exact=[s.sympify(z,locals=env)for z in record['exactCoefficients']]
  assert len(exact)==n+1 and all(z.free_symbols<=set(variables)for z in exact)
  vectors.append([s.simplify(z/s.sqrt(comb(n,r)))for r,z in enumerate(exact)])
  for z,saved in zip(exact,record['coefficients']):
   value=z.subs(subs);pair=[float(s.re(value)),float(s.im(value))]
   render_error=max(render_error,*(abs(a-b)for a,b in zip(pair,saved)))
 assert render_error<1e-15
 pivots=echelon(equations);assert len(pivots)==D
 components=0
 for a,V in enumerate(vectors):
  for b,W in enumerate(vectors):
   for i in range(3):
    for j in range(3):
     value=sum(comb(n-2,l)*V[i+l]*s.conjugate(W[j+l])for l in range(n-1))
     if a==b and i==j:value-=s.Rational(1,3*comb(2,i))
     for part in [s.re(value),s.im(value)]:
      assert in_span(poly(part,variables),pivots),(n,a,b,i,j)
      components+=1
 gram_components=0
 for a,V in enumerate(vectors):
  for b,W in enumerate(vectors):
   value=sum(comb(n,r)*V[r]*s.conjugate(W[r])for r in range(n+1))-int(a==b)
   for part in [s.re(value),s.im(value)]:assert in_span(poly(part,variables),pivots);gram_components+=1
 # Recompute derivative and Hessian bounds with independent rational arithmetic.
 residual=[evaluate(p,center)for p in equations]
 J=[[evaluate(derivative(p,j),center)for j in range(D)]for p in equations]
 H=[]
 for p in equations:
  h=[]
  for i in range(D):
   row=[]
   for j in range(D):
    dd=derivative(derivative(p,i),j)
    assert all(sum(powers)==0 for powers in dd)
    row.append(sum(dd.values()))
   h.append(row)
  H.append(h)
 eta=max(abs(dot(row,residual))for row in R)
 defect=[[Q(int(i==j))-sum(R[i][k]*J[k][j]for k in range(D))for j in range(D)]for i in range(D)]
 delta=maxrows(defect);Rnorm=maxrows(R);Hnorm=max(sum(abs(x)for row in h for x in row)for h in H)
 # This deliberately differs from N15's sharper transformed-Hessian bound.
 contraction=delta+radius*Rnorm*Hnorm
 assert delta<1 and contraction<1 and eta+contraction*radius<radius
 return {'n':n,'dimension':3,'orderAtLeast':2,'passed':True,'basisSha256':hashlib.sha256(basis_raw).hexdigest(),'certificateSha256':hashlib.sha256(cert_raw).hexdigest(),'variablesAndIndependentEquations':D,'allRealImaginaryCrossComponentsVerified':components,'separateGramComponentsVerified':gram_components,'previewCoefficientError':render_error,'boxRadius':str(radius),'newtonDisplacementStrictUpperBound':decimal_upper(eta),'jacobianDefectStrictUpperBound':decimal_upper(delta),'independentContractionStrictUpperBound':decimal_upper(contraction),'strictSelfMapping':True,'uniqueRealRootInClosedBox':True,'nonsingularJacobianThroughoutBox':True,'isolatedRootIsAlgebraic':True}

def verify():
 records=[audit(n)for n in[15,17]]
 return {'passed':True,'method':'Independent parsing of saved exact Dicke coefficients; rational row reduction checks every complete cross matrix against the defining equations. Fraction arithmetic independently reconstructs the Jacobian/Hessians and verifies a Newton-map Banach certificate with the separate bound ||I-RJ0||∞+r||R||∞ max_l sum_ij |H_lij|. No case-specific verifier or definition module is imported.','mathematicalJustification':['The closed rational box is complete. The Newton map x−RF(x) maps it strictly into itself and contracts, so it has exactly one fixed point.','The inverse defect is below one, making R and every Jacobian in the box nonsingular; a fixed point therefore satisfies F(x)=0.','A nonsingular isolated zero of a square rational polynomial system has algebraic coordinates.','Every full cross reduction equals delta_ab I3/3 at that root. Tracing gives an orthonormal basis, so every normalized vector in its span is second-order AC.'],'cases':records}
if __name__=='__main__':
 report=verify();(ROOT/'independent-certificate-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
