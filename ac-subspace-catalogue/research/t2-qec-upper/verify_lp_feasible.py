"""Independent standard-library rational audit of LP-feasible enumerators.
Uses factorial six-j sums, not the constructor's hypergeometric recurrence.
These certificates concern a relaxation, not existence of an AC2 subspace.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def transform(n):
 fac=[factorial(i)for i in range(2*n+2)]
 def entry(b,a):
  delta=F(fac[n-b]*fac[b]**2,fac[n+b+1])*F(fac[n-a]*fac[a]**2,fac[n+a+1])
  val=sum((F((-1)**z*fac[z+1],fac[z-n-b]**2*fac[z-n-a]**2*fac[2*n-z]*fac[n+a+b-z]**2)for z in range(n+max(a,b),min(2*n,n+a+b)+1)),F(0))
  return (-1)**(n+a+b)*(2*b+1)*delta*val
 return [[entry(b,a)for a in range(n+1)]for b in range(n+1)]

def verify(path):
 data=json.loads(path.read_text());n=data['n'];k=data['k'];d=n+1;x=list(map(F,data['point']));M=transform(n)
 dot=lambda u,v:sum((a*b for a,b in zip(u,v)),F(0))
 assert len(x)==d and sum(x)==1 and x[0]==F(k,d) and x[1]==x[2]==0
 assert dot(M[1],x)==dot(M[2],x)==0
 assert all(v>=0 for v in x)
 if n%2 and k%2:assert sum((-1)**a*x[a]for a in range(d))<=F(k-1,k)
 for l in range(d):
  B=dot(M[l],x);S=dot([(-1)**a*M[l][a]for a in range(d)],x)
  assert B>=0 and k*B>=x[l]
  assert 0<=S<=F(min(2*l+1,F(k*(k+(-1)**(n-l)),2)),k)
 c=F(n*(n+2),12);r=d-4*k
 fP=[F(1)]*d;fE=[3-F(a*(a+1),2)/c for a in range(d)];fR=[-4+F(a*(a+1),2)/c for a in range(d)]
 fE[0]=F(3);fR[0]=F(r,k);factors=[fP,fE,fR];dims=[k,3*k,r]
 for a in range(3):
  for b in range(a,3):
   weights=[x[l]*factors[a][l]*factors[b][l]for l in range(d)]
   assert sum(weights)==(F(dims[a],k)if a==b else 0)
   for l in range(d):
    B=dot(M[l],weights);S=dot([(-1)**h*M[l][h]for h in range(d)],weights)
    assert B>=0 and 0<=S<=F(2*l+1,k)
    if a==b:assert dims[a]*B>=weights[l]
 return {'n':n,'k':k,'passed':True,'pointFile':path.name,'independentFactorialSixJReconstruction':True,'allInequalitiesExact':True,'scope':'These necessary conditions cannot exclude the candidate dimension; no physical subspace is constructed.'}
if __name__=='__main__':
 result=[verify(HERE/f'lp-feasible-N{n}-k{k}.json')for n,k in [(35,8),(41,9),(42,9)]]
 (HERE/'lp-feasibility-independent-audit.json').write_text(json.dumps({'passed':True,'records':result},indent=2)+'\n');print(json.dumps(result,indent=2))
