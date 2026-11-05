#!/usr/bin/env python3
"""Independently verify all physical cross reductions of the saved coefficients."""
from pathlib import Path
import json, hashlib
import sympy as s
HERE=Path(__file__).resolve().parent

def verify(directory=None):
 folder=Path(directory) if directory is not None else HERE
 raw=(folder/'basis.json').read_bytes();d=json.loads(raw);n=d['n'];assert (n,d['t'],d['dimension'])==(19,3,2)
 vectors=[s.Matrix([s.sympify(x) for x in row['exactCoefficients']]) for row in d['basis']]
 assert all(len(c)==n+1 and not set().union(*(x.free_symbols for x in c)) for c in vectors)
 gram=s.Matrix.hstack(*vectors).H*s.Matrix.hstack(*vectors)
 assert all(s.simplify(s.expand_complex(x))==0 for x in gram-s.eye(2))
 M=[s.Matrix(4,n-2,lambda a,b:c[a+b]*s.sqrt(s.binomial(3,a)*s.binomial(n-3,b)/s.binomial(n,a+b))) for c in vectors]
 count=0
 for i,A in enumerate(M):
  for j,B in enumerate(M):
   for value in A*B.H-int(i==j)*s.eye(4)/4:
    assert s.simplify(s.expand_complex(value))==0;count+=1
 error=max(abs(complex(s.N(z,30))-complex(*approx)) for c,row in zip(vectors,d['basis']) for z,approx in zip(c,row['coefficients']));assert error<1e-12
 return {'passed':True,'n':n,'t':3,'dimension':2,'orthonormalityVerified':True,'allOrderedCrossReductionsVerified':True,'complexCrossEntriesExactlyVerified':count,'fullCrossReductionEntriesChecked':count,'basisSHA256':hashlib.sha256(raw).hexdigest(),'maximumDisplayCoefficientError':error,'method':'Read exported radical coefficients, construct the normalized-Dicke 3|(N-3) coefficient matrices, and symbolically verify every entry of all four ordered reduced matrices equals delta_ab I4/4. No constructor or primary verifier is imported.'}
if __name__=='__main__':
 result=verify();(HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
