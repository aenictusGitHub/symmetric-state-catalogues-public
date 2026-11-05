#!/usr/bin/env python3
"""Independent physical-spin verification of an exported radical AC basis.
No primary verifier, constructor or numerical source is imported.
"""
from pathlib import Path
import json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent
EXPECTED=(18, 2, 4)

def verify():
 raw=(HERE/'basis.json').read_bytes();b=json.loads(raw)
 n,t,k=b['n'],b['t'],b['dimension'];assert (n,t,k)==EXPECTED
 B=s.Matrix.hstack(*[s.Matrix([s.sympify(v) for v in row['exactCoefficients']]) for row in b['basis']])
 assert B.shape==(n+1,k) and all(not v.free_symbols for v in B)
 Z=s.diag(*[s.Rational(n,2)-r for r in range(n+1)]);J=s.zeros(n+1)
 for r in range(1,n+1):J[r-1,r]=s.sqrt(r*(n-r+1))
 count=0;operators=[]
 for q in range(t+1):
  for p in range(t-q+1):
   A=J**q*Z**p;expected=s.trace(A)/s.Integer(n+1)
   C=B.H*A*B-expected*s.eye(k)
   for value in C:
    assert s.simplify(s.expand_complex(value))==0,(q,p,value);count+=1
   operators.append({'raisingPower':q,'zPower':p,'scalarExpectation':str(expected)})
 maximum=0
 for row in b['basis']:
  for v,approx in zip(row['exactCoefficients'],row['coefficients']):
   error=abs(complex(s.N(s.sympify(v),30))-complex(*approx));maximum=max(maximum,error);assert error<3e-12
 return {'passed':True,'n':n,'t':t,'dimension':k,'maximumDimension':k,'orthonormalityVerified':True,'physicalComplexMomentEntriesChecked':count,'operators':operators,'complexCrossEntriesExactlyVerified':k*k*(t+1)**2,'fullCrossReductionEntriesChecked':k*k*(t+1)**2,'basisSHA256':hashlib.sha256(raw).hexdigest(),'maximumDisplayCoefficientError':maximum,'method':'Independently construct physical Jz and J+ matrices and check every ordered basis-pair compression of J+^q Jz^p for q+p<=t. These operators and their adjoints span all spin-tensor ranks zero through t; their scalar compressions are equivalent to every ordered t-qubit cross reduction. No primary verifier or constructor is imported.','coefficientType':'Explicit finite radical expressions, with no free parameters.'}

if __name__=='__main__':
 result=verify();(HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
