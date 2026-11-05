#!/usr/bin/env python3
"""Exact N18 AC5 witness and dimension bound; Python standard library only."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
from collections import defaultdict
import hashlib,json
HERE=Path(__file__).resolve().parent

def sqrt_fraction(q):
 q=F(q);assert q>=0
 if not q:return F(0),1
 v=q.numerator*q.denominator;square=1;radicand=1;p=2
 while p*p<=v:
  count=0
  while v%p==0:v//=p;count+=1
  square*=p**(count//2)
  if count%2:radicand*=p
  p+=1
 if v>1:radicand*=v
 return F(square,q.denominator),radicand

def reduction(t):
 n=18;weights={0:F(14,130),5:F(51,130),13:F(51,130),18:F(14,130)};phases={0:0,5:1,13:1,18:0};matrix=[]
 for a in range(t+1):
  row=[]
  for b in range(t+1):
   entry=defaultdict(F)
   for l in range(n-t+1):
    r,s=a+l,b+l
    if r not in weights or s not in weights:continue
    squared=weights[r]*weights[s]*F(comb(t,a)*comb(t,b)*comb(n-t,l)**2,comb(n,r)*comb(n,s))
    coefficient,radicand=sqrt_fraction(squared);phase=(phases[r]-phases[s])%4
    entry[phase%2,radicand]+=coefficient*(1 if phase<2 else -1)
   row.append({key:value for key,value in entry.items() if value})
  matrix.append(row)
 return matrix

def distance_squared(matrix):
 dim=len(matrix);total=defaultdict(F)
 for a,row in enumerate(matrix):
  for b,entry in enumerate(row):
   entry=dict(entry)
   if a==b:entry[0,1]=entry.get((0,1),F(0))-F(1,dim)
   for (component,rad),value in entry.items():
    for (component2,rad2),value2 in entry.items():
     if component!=component2:continue
     factor,radicand=sqrt_fraction(rad*rad2);total[radicand]+=value*value2*factor
 assert all(not v for r,v in total.items() if r!=1)
 return total.get(1,F(0))

def verify():
 weights={0:F(14,130),5:F(51,130),13:F(51,130),18:F(14,130)}
 assert sum(weights.values())==1
 # Source entries are 2sqrt455, i sqrt6630, i sqrt6630, 2sqrt455.
 assert 2*(4*455+6630)==130**2
 assert F(4*455,130**2)==weights[0] and F(6630,130**2)==weights[5]
 moments=[]
 for p in range(6):
  value=sum(w*(i-9)**p for i,w in weights.items());target=F(sum((i-9)**p for i in range(19)),19);assert value==target
  moments.append({'power':p,'value':str(value)})
 reductions=[]
 for t in range(1,6):
  rho=reduction(t)
  for a,row in enumerate(rho):
   for b,entry in enumerate(row):assert entry==({(0,1):F(1,t+1)}if a==b else {})
  assert distance_squared(rho)==0
  reductions.append({'order':t,'entriesChecked':(t+1)**2,'identityVerified':True})
 sixth=distance_squared(reduction(6));assert sixth==F(73,59150)>0
 # Complete rectangular irreps from spin9 tospin17/2: ranks1/2,3/2,5/2.
 ranks=[F(1,2),F(3,2),F(5,2)];count=sum(int(2*l+1)for l in ranks);assert count==12
 upper=18//count;assert upper==1
 report={'passed':True,'n':18,'t':5,'maximumDimension':1,'stateLatex':r'|\psi\rangle=\frac{\sqrt{14}(|D_{18}^{(0)}\rangle+|D_{18}^{(18)}\rangle)+i\sqrt{51}(|D_{18}^{(5)}\rangle+|D_{18}^{(13)}\rangle)}{\sqrt{130}}',
 'source':{'file':'ACLiege.jl','line':203,'attribution':'M. Grassl, as credited by the source comment','sameSourceRay':True},
 'probabilities':{str(i):str(w)for i,w in weights.items()},'exactDiagonalMoments':moments,'allReducedMatrices':reductions,'exactReducedMatrixEntriesChecked':sum((t+1)**2 for t in range(1,6)),
 'sixthOrderFrobeniusDistanceSquared':str(sixth),'exactAchievedOrder':5,
 'upperBound':{'rectangularRanks':[str(l)for l in ranks],'imagesPerBasisVector':count,'outputDimension':18,'inequality':'12k<=18','wholeSubspacePolarizationRequired':True,'upper':upper},
 'method':'Independent normalized-Dicke splitting with exact fractions and squarefree radicals; every reduced matrix entry through order5, plus the next-order distance. Standard library only.',
 'scope':'Exact witness and optimal fifth-order subspace dimension; the witness calculation alone does not exclude a different sixth-order state.'}
 return report

if __name__=='__main__':
 result=verify();(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
