#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent

def verify():
 u=s.zeros(14,1);v=s.zeros(14,1);u[0]=s.sqrt(10)/6;u[9]=s.sqrt(26)/6;v[4]=-s.sqrt(26)/6;v[13]=s.sqrt(10)/6
 B=s.Matrix.hstack(u,v);assert (B.H*B-s.eye(2)).applyfunc(s.simplify)==s.zeros(2)
 def split(c,t):return s.Matrix(t+1,14-t,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(13-t,b)/s.binomial(13,a+b)))
 M=[split(c,2)for c in[u,v]];count=0
 for a,A in enumerate(M):
  for b,D in enumerate(M):
   for z in A*D.H-s.eye(3)*s.Rational(int(a==b),3):assert s.simplify(z)==0;count+=1
 counter=(split(u,3)*split(u,3).H-s.eye(4)/4).applyfunc(s.simplify);assert counter!=s.zeros(4)
 return {'passed':True,'n':13,'dimension':2,'order':2,'wholeSubspaceCrossEntriesChecked':count,'orthonormal':True,'exactCoefficients':[[str(z)for z in c]for c in[u,v]],'thirdOrderCounterexample00':str(counter[0,0])}
if __name__=='__main__':
 r=verify();(ROOT/'lower-verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
