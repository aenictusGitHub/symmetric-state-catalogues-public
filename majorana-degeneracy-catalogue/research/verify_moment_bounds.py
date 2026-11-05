#!/usr/bin/env python3
"""Verify moment feasibility/exclusion certificates with rational arithmetic."""
from pathlib import Path
from fractions import Fraction as F
import json

def evaluate(p,k):
 v=F(0)
 for c in reversed(p):v=v*k+c
 return v

def main():
 path=Path(__file__).with_name('moment-bounds.json');d=json.loads(path.read_text());counts={'moment_excluded':0,'moment_feasible':0}
 for r in d['records']:
  n,g,t=r['n'],r['g'],r['t'];counts[r['status']]+=1
  if r['status']=='moment_excluded':
   p=list(map(F,r['polynomialAscending']));assert len(p)-1<=t
   assert all(evaluate(p,k)>=0 for k in range(g,n+1))
   expected=sum(evaluate(p,k) for k in range(n+1))/F(n+1)
   assert expected==F(r['uniformExpectation'])<0
  else:
   support=r['support'];weights=list(map(F,r['weights']));assert len(support)==len(weights)
   assert all(g<=k<=n for k in support) and all(w>=0 for w in weights)
   for degree in range(t+1):
    assert sum(w*k**degree for k,w in zip(support,weights))==sum(F(k**degree,n+1) for k in range(n+1))
 print(json.dumps({'records':len(d['records']),'counts':counts,'exactRationalVerification':True},indent=2))

if __name__=='__main__':main()
