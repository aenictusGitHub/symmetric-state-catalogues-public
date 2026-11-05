#!/usr/bin/env python3
"""Verify two exact12-qubit full-partition witnesses; requires SymPy only."""
import json
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;z=s.Symbol('z');plus=z**4+2*s.I*s.sqrt(3)*z**2+1;minus=z**4-2*s.I*s.sqrt(3)*z**2+1
A=s.Poly(plus,z,extension=True);B=s.Poly(minus,z,domain=A.domain)
assert s.gcd(A,A.diff()).degree()==s.gcd(B,B.diff()).degree()==s.gcd(A,B).degree()==0
expected=[(plus**3,[3,3,3,3],s.Rational(441,14884)),(plus**2*minus,[2,2,2,2,1,1,1,1],s.Rational(1,484))];records=json.loads((HERE/'states.json').read_text())['states'];assert len(records)==2
results=[]
for data,(polynomial,partition,distance) in zip(records,expected):
 c=list(map(s.sympify,data['exactCoefficients']));assert data['n']==12 and len(c)==13
 assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
 f=s.Poly(sum(s.simplify((-1)**r*s.sqrt(s.binomial(12,r))*x)*z**r for r,x in enumerate(c)),z,extension=True)
 assert s.simplify(f.monic().as_expr()-s.expand(polynomial))==0
 assert data['partition']==partition and sum(partition)==12
 for t in (1,2,3):
  M=s.zeros(t+1,13-t)
  for a in range(t+1):
   for q in range(13-t):M[a,q]=c[a+q]*s.sqrt(s.binomial(t,a)*s.binomial(12-t,q)/s.binomial(12,a+q))
  rho=(M*s.conjugate(M.T)).applyfunc(s.simplify);delta=rho-s.eye(t+1)/s.Integer(t+1);sq=s.simplify(s.trace(delta*delta))
  assert sq==(distance if t==3 else 0)
 results.append({'id':data['id'],'partition':partition,'exactWitnessOrder':2,'rho2IdentityVerified':True,'rho3DistanceSquared':str(distance),'partitionLowerBound':2,'partitionMaximumNotInferredFromWitness':True})
print(json.dumps({'verified':True,'method':'Exact normalized Dicke coefficients, reduced matrices and coprime squarefree tetrahedral quartics','states':results},indent=2))
