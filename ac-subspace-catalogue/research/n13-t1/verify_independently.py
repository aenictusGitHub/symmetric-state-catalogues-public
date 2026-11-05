#!/usr/bin/env python3
"""Independent spin-word verification of the N=13,k=6,t=1 exclusion.
Uses exact integer/Fraction arithmetic. No optimization, Wigner symbols or
constructor imports are used; data records provide basis/certificate only.
"""
from pathlib import Path
from itertools import combinations
from math import comb, factorial, gcd,lcm
from fractions import Fraction as F
import json,hashlib,time,gzip

N=13;K=6

def canonical(word):
 if len(set(word))<len(word):return None,0
 inv=sum(word[a]>word[b] for a in range(len(word))for b in range(a+1,len(word)))
 return tuple(sorted(word)),(-1)**inv

def product(xs):
 r=1
 for x in xs:r*=x
 return r

def gram(U,ws):
 m=len(U[0]);return [[sum(product(comb(N,r)for r in w)*U[i][a]*U[i][b]for i,w in enumerate(ws))for b in range(m)]for a in range(m)]

def move(U,ws,direction,dst):
 m=len(U[0]);out=[[0]*m for _ in dst];loc={w:i for i,w in enumerate(dst)}
 for a,w in enumerate(ws):
  for leg,r in enumerate(w):
   rr=r+direction
   if not 0<=rr<=N:continue
   z=list(w);z[leg]=rr;q,sg=canonical(z)
   if q is None:continue
   factor=(r+1)if direction==1 else N-r+1
   target=out[loc[q]]
   for b in range(m):target[b]+=sg*factor*U[a][b]
 return out

def rankmod(A,p=1000000007):
 A=[[int(x)%p for x in row]for row in A];r=0
 for c in range(len(A[0])):
  pivot=next((a for a in range(r,len(A))if A[a][c]),None)
  if pivot is None:continue
  A[r],A[pivot]=A[pivot],A[r];inv=pow(A[r][c],p-2,p);A[r]=[(x*inv)%p for x in A[r]]
  for a in range(r+1,len(A)):
   z=A[a][c]
   if z:A[a]=[(x-z*y)%p for x,y in zip(A[a],A[r])]
  r+=1
  if r==len(A):break
 return r

def primitive(row):
 den=lcm(*(x.denominator for x in row));out=[int(x*den)for x in row];g=gcd(*out)
 return [x//g for x in out]if g else out

def verify_basis(blocks):
 by={}
 for w in combinations(range(N+1),K):by.setdefault(2*sum(w)-K*N,[]).append(w)
 descendants=[];gram_checks=0
 for b in blocks:
  J=b['J2'];H=b['highest'];m=b['multiplicity']
  assert m==len(by[J])-len(by.get(J+2,[]))
  assert len(H)==len(by[J])and all(len(row)==m for row in H)
  assert rankmod(H)==m
  raised=move(H,by[J],1,by.get(J+2,[]));assert not any(x for row in raised for x in row)
  G=gram(H,by[J]);U={J:H}
  if 'norms'in b:assert b['norms']==[factorial(K)*G[a][a]for a in range(m)]
  for step in range(1,J+1):
   M=J-2*step;raw=move(U[M+2],by[M+2],-1,by[M]);assert all(x%step==0 for row in raw for x in row)
   U[M]=[[x//step for x in row]for row in raw];GM=gram(U[M],by[M])
   assert GM==[[comb(J,step)*x for x in row]for row in G];gram_checks+=m*m
  descendants.append(U)
 assert sum((b['J2']+1)*b['multiplicity']for b in blocks)==comb(N+1,K)==3003
 # Exact weight-character equality ensures the decomposition is complete.
 for M,ws in by.items():assert len(ws)==sum(b['multiplicity']for b in blocks if b['J2']>=abs(M))
 return by,descendants,gram_checks

def physical_row(meta,blocks,variables,by,Us):
 # Metadata supports [weight_of_fourword,alpha_index,beta_index].
 Ml,ai,bi=map(int,meta)
 lesser=[w for w in combinations(range(N+1),K-1)if 2*sum(w)-(K-1)*N==Ml]
 alpha,beta=lesser[ai],lesser[bi];idx={M:{w:i for i,w in enumerate(ws)}for M,ws in by.items()}
 row=[]
 for block,a,b in variables:
  J=blocks[block]['J2'];z=F(0)
  for i in range(N+1):
   M=Ml+2*i-N
   if M not in Us[block]:continue
   wa,sa=canonical((i,)+alpha);wb,sb=canonical((i,)+beta)
   if wa is None or wb is None:continue
   ua=Us[block][M][idx[M][wa]];ub=Us[block][M][idx[M][wb]]
   value=ua[a]*ub[b]+(ua[b]*ub[a]if a!=b else 0)
   z+=F(sa*sb*comb(N,i)*(2*i-N)*value,comb(J,(J-M)//2))
  row.append(z)
 return primitive(row)

def verify(directory=None):
 start=time.time();root=Path(directory)if directory is not None else Path(__file__).resolve().parent
 inputs={name:((root/name).read_bytes() if (root/name).exists() else gzip.decompress((root/(name+'.gz')).read_bytes())) for name in ['kernel.json','certificate.json']}
 data=json.loads(inputs['kernel.json']);certificate=json.loads(inputs['certificate.json'])
 assert(data['n'],data['k'],data['order'],data['exteriorPower'],data['dimension'])==(13,6,1,6,3003)
 bs=data['blocks'];variables=data['variables'];by,Us,gram_checks=verify_basis(bs)
 assert variables==[[ib,a,b]for ib,B in enumerate(bs)for a in range(B['multiplicity'])for b in range(a,B['multiplicity'])]
 assert len(variables)==735
 rebuilt=[physical_row(meta,bs,variables,by,Us)for meta in data['selectedConstraints']]
 assert rebuilt==data['selectedRows']
 assert len(rebuilt)==716 and rankmod(rebuilt)==716
 kernel=data['kernel'];assert len(kernel)==735 and all(len(row)==19 for row in kernel)
 assert rankmod(kernel)==19
 assert all(sum(row[i]*kernel[i][j]for i in range(735))==0 for row in rebuilt for j in range(19))
 diagonals=[i for i,(_,a,b)in enumerate(variables)if a==b]
 assert certificate['diagonalVariableIndices']==diagonals and len(diagonals)==151
 weights=certificate['positiveCoefficients'];assert len(weights)==151 and all(isinstance(x,int)and x>0 for x in weights)
 assert all(sum(weights[a]*kernel[index][j]for a,index in enumerate(diagonals))==0 for j in range(19))
 report={'passed':True,'n':13,'order':1,'excludedDimension':6,'method':'Independent exact signed exterior-word contractions, with no constructor imports, optimizer or Wigner functions.',
 'inputSha256':{name:hashlib.sha256(raw).hexdigest()for name,raw in inputs.items()},
 'exteriorDimension':3003,'highestWeights':151,'spinBlocks':24,'realSymmetricInvariantCoordinates':735,
 'descendantGramEntriesChecked':gram_checks,'reconstructedPhysicalRows':716,'exactMatrixKernelProductsChecked':13604,
 'constraintRank':716,'kernelRank':19,'rankPrime':1000000007,
 'strictlyPositiveDiagonalCoefficients':151,'exactPositiveFunctionalKernelProductsChecked':19,
 'largestCoefficientDecimalDigits':max(len(str(x))for x in weights),
 'complexScope':'For any complex code, the real part of every positive invariant multiplicity block is positive and obeys these real contraction equations. Positive diagonal cancellation forces all blocks, including imaginary off-diagonal entries, to vanish.',
 'conclusion':'No six-dimensional first-order anticoherent subspace exists in spin 13/2.',
 'elapsedSeconds':round(time.time()-start,3)}
 (root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
 return report

if __name__=='__main__':
 import sys
 print(json.dumps(verify(sys.argv[1]if len(sys.argv)>1 else None),indent=2))
