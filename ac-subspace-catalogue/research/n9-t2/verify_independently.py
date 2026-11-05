#!/usr/bin/env python3
"""Independent exact physical-index audit of the N9,t2,k2 Young22 obstruction."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
from fractions import Fraction as F
from math import comb,gcd,lcm,prod
import json,time,hashlib
import sympy as s

def verify(directory=None):
 ROOT=Path(directory) if directory is not None else Path(__file__).resolve().parent;SRC=ROOT/'kernel.json';raw=SRC.read_bytes();data=json.loads(raw);n=data['n'];d=n+1;start=time.time()
 pairs=list(combinations(range(d),2));words=[p+q for p in pairs for q in pairs if p<=q]
 by=defaultdict(list)
 for w in words:by[2*sum(w)-4*n].append(w)
 weights={w:prod(comb(n,r) for r in w) for w in __import__('itertools').product(range(d),repeat=4)}
 def orbit(w):
  a,b,c,e=w;out={}
  for left,right in [((a,b),(c,e)),((c,e),(a,b))]:
   for p,sp in [(left,1),(left[::-1],-1)]:
    for q,sq in [(right,1),(right[::-1],-1)]:
     word=p+q;sign=sp*sq
     assert word not in out or out[word]==sign
     out[word]=sign
  return out

 def step_physical(vector,direction,divisor=1):
  out=defaultdict(int)
  for source,value in vector.items():
   for leg,r in enumerate(source):
    target=list(source);target[leg]+=direction
    if not 0<=target[leg]<=n:continue
    # Acting directly on all physical amplitudes, with sqrt(binomial) weights
    # removed. Direction-1 is physical J_+, direction+1 is physical J_-.
    factor=n-r+1 if direction==-1 else r+1
    out[tuple(target)]+=factor*value
  result={}
  for word,val in out.items():
   assert val%divisor==0
   if val:result[word]=val//divisor
  return result

 def inner(u,v):return sum(value*v.get(w,0)*weights[w] for w,value in u.items())
 # Independently count the SU(2) character of Sym^2(Lambda^2 V)-Lambda^4 V.
 # This proves the supplied highest-weight list is complete, not just independent.
 quadweights=Counter(2*sum(q)-4*n for q in combinations(range(d),4))
 weightdim={M:len(by.get(M,[]))-quadweights[M] for M in range(0,4*n+3,2)}
 expected={M:weightdim[M]-weightdim.get(M+2,0) for M in range(0,4*n+1,2) if weightdim[M]-weightdim.get(M+2,0)}
 actual={b['J2']:b['multiplicity'] for b in data['blocks']}
 assert len(actual)==len(data['blocks']) and actual==expected
 assert sum((M+1)*mult for M,mult in expected.items())==d*d*(d*d-1)//12==825
 blocks=[];highest_checks=0;norm_checks=0
 for ib,block in enumerate(data['blocks']):
  M=block['J2'];mult=block['multiplicity'];highest=block['highest'];physical=[]
  assert len(highest)==len(by[M])
  assert s.Matrix(highest).rank()==mult
  for a in range(mult):
   vec={}
   for word,row in zip(by[M],highest):
    if row[a]:
     for w,sign in orbit(word).items():
      assert w not in vec;vec[w]=sign*row[a]
   assert not step_physical(vec,1)
   for q in combinations(range(d),4):
    if 2*sum(q)-4*n==M:
     a0,b0,c0,e0=q
     assert vec.get((a0,b0,c0,e0),0)-vec.get((a0,c0,b0,e0),0)+vec.get((a0,e0,b0,c0),0)==0
   assert inner(vec,vec)==block['norms'][a]
   physical.append(vec);highest_checks+=1
  gram=[[inner(a,b) for b in physical] for a in physical]
  U={M:physical}
  for k in range(1,M+1):
   physical=[step_physical(v,-1,k) for v in physical];U[M-2*k]=physical
   for a in range(mult):
    for b in range(a,mult):
     assert inner(physical[a],physical[b])==comb(M,k)*gram[a][b];norm_checks+=1
  blocks.append({'J2':M,'U':U,'gram':gram})
 assert sum((b['J2']+1)*data['blocks'][i]['multiplicity'] for i,b in enumerate(blocks))==825
 low=defaultdict(list)
 for a in range(d):
  for b,c in pairs:low[2*(a+b+c)-3*n].append((a,b,c))
 A=[];row_scales=[]
 for rowindex,(ell,M3,alpha,beta) in enumerate(data['selectedConstraints']):
  left=low[M3][alpha];right=low[M3][beta];row=[]
  for ib,a,b in data['variables']:
   block=blocks[ib];J2=block['J2'];value=F(0)
   for i in range(d):
    M=M3+2*i-n
    if M not in block['U']:continue
    v=block['U'][M];f=(2*i-n) if ell==1 else 3*(2*i-n)**2-n*(n+2)
    product=v[a].get((i,)+left,0)*v[b].get((i,)+right,0)
    if a!=b:product+=v[b].get((i,)+left,0)*v[a].get((i,)+right,0)
    value+=F(comb(n,i)*f*product,comb(J2,(J2-M)//2))
   row.append(value)
  den=lcm(*(x.denominator for x in row));primitive=[x.numerator*(den//x.denominator) for x in row];common=gcd(*primitive);assert common
  primitive=[x//common for x in primitive];row_scales.append(str(F(den,common)));A.append(primitive)
 # Match discovery rows up to their harmless nonzero primitive scales.
 assert data['variables']==[[ib,a,b] for ib,block in enumerate(data['blocks']) for a in range(block['multiplicity']) for b in range(a,block['multiplicity'])]
 assert len(data['selectedRows'])==156
 for physical,saved in zip(A,data['selectedRows']):
  pivot=next(i for i,x in enumerate(physical) if x)
  assert saved[pivot]!=0
  assert all(x*saved[pivot]==y*physical[pivot] for x,y in zip(physical,saved))
 K=[[int(v) for v in row] for row in data['kernel']]
 assert len(K)==169 and all(len(row)==13 for row in K)
 for row in A:
  for col in range(13):assert sum(x*K[i][col] for i,x in enumerate(row))==0

 def rank_mod(rows,prime=1000003):
  assert s.isprime(prime)
  work=[[int(v)%prime for v in row] for row in rows];rank=0
  for col in range(len(work[0])):
   pivot=next((r for r in range(rank,len(work)) if work[r][col]),None)
   if pivot is None:continue
   work[rank],work[pivot]=work[pivot],work[rank]
   inverse=pow(work[rank][col],-1,prime)
   work[rank]=[(v*inverse)%prime for v in work[rank]]
   for r in range(rank+1,len(work)):
    if work[r][col]:
     scale=work[r][col];work[r]=[(v-scale*u)%prime for v,u in zip(work[r],work[rank])]
   rank+=1
   if rank==len(work):break
  return rank
 assert rank_mod(A)==156 and rank_mod(K)==13
 cert=json.loads((ROOT/'certificate.json').read_text())
 diagonal=[i for i,(ib,a,b) in enumerate(data['variables']) if a==b]
 assert len(diagonal)==57 and diagonal==cert['diagonalVariableIndices']
 positive=[int(v) for v in cert['positiveCoefficients']]
 assert len(positive)==57 and all(v>0 for v in positive)
 for col in range(13):assert sum(c*K[i][col] for i,c in zip(diagonal,positive))==0
 report={'passed':True,'n':9,'t':2,'excludedDimension':2,'sourceSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent full physical-word expansion, physical spin-ladder operations and direct first-leg contractions; no canonical ladder or contraction code imported. Exact integer equations and modular full-rank certificates prove the positive functional vanishes on the complete constraint kernel.','highestWeightVectorsIndependentlyChecked':highest_checks,'loweringCrossGramIdentitiesIndependentlyChecked':norm_checks,'young22Dimension':825,'schurMultiplicitiesIndependentlyCounted':True,'constraintRowsReconstructed':156,'variables':169,'exactConstraintRank':156,'exactKernelRank':13,'strictlyPositiveDiagonalCoefficients':57,'allKernelProductsExactlyZero':True,'allSelectedPhysicalRowsAgreeWithSource':True,'noFloatingPointProofSteps':True,'realSymmetricRestriction':'Without loss: average the SU(2)-invariant moment with its entrywise conjugate. The Young support, all rank1/rank2 constraints, trace and positivity are preserved. The real highest-weight basis turns this into entrywise conjugation of every invariant block.','positiveSemidefiniteArgument':'Every invariant coefficient block is positive semidefinite even though its multiplicity basis is not orthogonal: the full intertwiner is injective, so positivity is equivalent under congruence. All57 diagonal coordinates are nonnegative. Their strictly positive weighted sum vanishes, forcing every diagonal and hence every block to vanish. This contradicts trace1 of the twirled normalized wedge-square moment.'}
 return report

if __name__=='__main__':
 report=verify();(Path(__file__).with_name('independent-verification.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
