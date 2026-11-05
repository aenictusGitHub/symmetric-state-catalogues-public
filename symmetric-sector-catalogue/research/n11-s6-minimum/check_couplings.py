"""Verify every transvectant intertwiner using exact integer recurrences."""
from pathlib import Path
from math import comb,factorial
import json
R=Path(__file__).parent;d=json.load(open(R/'covariants.json'));raw=json.load(open(R/'source-basis.json'))['blocks'];seen=set()
def ff(a,b):return factorial(a)//factorial(a-b)if a>=b>=0 else 0
for B in d['blocks']:
 L=B['L2'];p,q=map(int,B['kind']);assert len(B['channels'])==len(B['mixing'])==B['dimension'];assert all(len(row)==B['dimension']for row in B['mixing'])
 for ch in B['channels']:
  ib=ch['sourceBlock'];J=raw[ib]['L2'];ps,qs=map(int,raw[ib]['kind']);assert ib in d['sourceBlocks'];assert (p,q)==((ps+1,qs)if ch['coupledWith']=='f'else(ps,qs+1));assert 0<=ch['sourceChannel'];assert abs(J-11)<=L<=J+11 and(J+11-L)%2==0
  if(J,L)in seen:continue
  seen.add((J,L));s=(J+11-L)//2;T={}
  for a in range(J+1):
   for b in range(12):
    r=a+b-s
    if 0<=r<=L:T[r,a,b]=sum((-1)**k*comb(s,k)*ff(J-a,s-k)*ff(a,k)*ff(11-b,k)*ff(b,s-k)for k in range(s+1))
  for r in range(L+1):
   for a in range(J+1):
    for b in range(12):
     assert a*T.get((r,a-1,b),0)+b*T.get((r,a,b-1),0)==(r+1)*T.get((r+1,a,b),0)
     assert (J-a)*T.get((r,a+1,b),0)+(11-b)*T.get((r,a,b+1),0)==(L-r+1)*T.get((r-1,a,b),0)
print('PASS exact raising/lowering identities for',len(seen),'transvectants')
