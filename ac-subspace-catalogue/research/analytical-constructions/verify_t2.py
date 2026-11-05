#!/usr/bin/env python3
"""Exact cross-reduction verification of the independent N=11 and N=16 constructions."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent

def split(c,n,t):
 return s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b)))
def verify(record):
 n=record['n'];t=record['t'];C=[s.Matrix([s.sympify(c)for c in b['exactCoefficients']])for b in record['basis']];B=s.Matrix.hstack(*C)
 assert s.simplify(B.H*B)==s.eye(len(C))
 M=[split(c,n,t)for c in C];count=0
 for a,A in enumerate(M):
  for b,D in enumerate(M):
   for z in A*D.H-s.eye(t+1)*s.Rational(int(a==b),t+1):
    assert s.simplify(s.expand_complex(z))==0
    count+=1
 nextM=[split(c,n,t+1)for c in C];next_bad=any(s.simplify(z)!=0 for a,A in enumerate(nextM)for b,D in enumerate(nextM)for z in A*D.H-s.eye(t+2)*s.Rational(int(a==b),t+2));assert next_bad
 out={'id':record['id'],'passed':True,'exactFullCrossEntries':count,'orthonormal':True,'achievedOrder':t,'nextOrderFails':True,'dimension':len(C)}
 if n==11:
  j=s.Rational(n,2);factors=[s.factor(4-3*L*(L+1)/(2*j*(j+1)))for L in range(1,n+1)]
  assert all(factors)
  out.update(maximumDimension=2,upperBound='Tensor-Gram saturation: k=3 fills all 12 dimensions, forcing a projector P=I/4 because every non-scalar tensor factor is nonzero. This is impossible.',nonzeroTensorFactors=[str(v)for v in factors])
 if n==16:
  q=s.sympify(record['phase']);assert s.simplify(q*s.conjugate(q))==1
  out.update(rigorousDimensionRange=[3,4],upperBound='Tensor-Gram bound floor((N+1)/4). The witness alone does not exclude dimension four.')
 return out
if __name__=='__main__':
 data=json.loads((ROOT/'t2-bases.json').read_text())
 reports=[verify(r)for r in data['subspaces']]
 report={'passed':True,'method':'Exact Dicke splitting and every ordered cross-reduction block, including cross terms between distinct basis states.','witnesses':reports}
 (ROOT/'t2-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
