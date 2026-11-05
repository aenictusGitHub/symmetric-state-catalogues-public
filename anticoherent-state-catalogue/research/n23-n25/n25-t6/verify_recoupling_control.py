from fractions import Fraction as Q
from math import comb,prod
from functools import lru_cache
from collections import Counter
from pathlib import Path
import sys,json,time
from coupled_engine import cg_highest,sixj_raw
start=time.time()
@lru_cache(None)
def descend(n,v_tuple,step):
 v=dict(v_tuple);out={}
 for w,c in v.items():
  for p,i in enumerate(w):
   if i:
    z=list(w);z[p]-=1;z=tuple(z);out[z]=out.get(z,0)+c*(n-i+1)
 assert all(c%step==0 for c in out.values())
 return tuple((w,c//step)for w,c in out.items()if c)
@lru_cache(None)
def pair(n,P,s):
 k,lo,cs,_=cg_highest(n,n,P)
 v=tuple(((n-q,n-k+q),c)for q,c in zip(range(lo,lo+len(cs)),cs))
 for d in range(1,s+1):v=descend(n,v,d)
 return v
@lru_cache(None)
def coupled(n,P,S,J,step=0):
 k,lo,cs,_=cg_highest(P,S,J);out={}
 for q,c in zip(range(lo,lo+len(cs)),cs):
  for a,aa in pair(n,P,q):
   for b,bb in pair(n,S,k-q):out[a+b]=out.get(a+b,0)+c*aa*bb
 v=tuple((w,c)for w,c in out.items()if c)
 for d in range(1,step+1):v=descend(n,v,d)
 return v
@lru_cache(None)
def met(n,w):return prod(comb(n,i)for i in w)
results=[]
for n in [1,3,5]:
 ps=list(range(2*n,-1,-4));count=0
 for P in ps:
  for S in ps:
   for R in ps:
    for T in ps:
     norms=prod(cg_highest(n,n,p)[3]for p in[P,S,R,T])
     for J in range(max(abs(P-S),abs(R-T)),min(P+S,R+T)+1,2):
      for L in range(max(abs(R-S),abs(P-T)),min(R+S,P+T)+1,2):
       u=dict(coupled(n,R,S,L));v=dict(coupled(n,P,T,L));val=Q(0)
       for s in range(J+1):
        partial=0
        for a,aa in coupled(n,P,S,J,s):
         for b,bb in coupled(n,R,T,J,s):
          row=(n-b[0],n-b[1],a[2],a[3]);col=(n-a[0],n-a[1],b[2],b[3]);uv=u.get(row,0)*v.get(col,0)
          if uv:partial+=(-1)**sum(a[:2]+b[:2])*uv*aa*bb*met(n,row)*met(n,col)
        val+=Q(partial,comb(J,s))
       val*=L+1
       expected=(J+1)*(L+1)*sixj_raw(P,S,J,R,T,L)*norms
       assert val==expected,(n,P,S,R,T,J,L,val,expected)
       count+=1
 results.append(dict(n=n,exactCases=count));print(results[-1],flush=True)
Path(__file__).resolve().with_name('recoupling-control.json').write_text(json.dumps(dict(passed=True,tests=results,method='Independent ordered four-copy amplitudes in the binomial basis, exact Fraction arithmetic, explicit first-pair partial time reversal and physical metric.',seconds=time.time()-start),indent=2)+'\n')
