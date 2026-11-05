#!/usr/bin/env python3
"""Exact conditional N=13,t=2 exclusion. Python standard library only.
Reconstructs all physical spin actions and exterior contractions independently
of the optimization/constructor code. The certificate is checked, not trusted.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,gcd,lcm
import json,hashlib,time
N=13;D=14
# Exact real quadratic field Q(sqrt(11)), represented by pairs.
def f(a=0,b=0):return(F(a),F(b))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def mul(x,y):return(x[0]*y[0]+11*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,a):return(x[0]*a,x[1]*a)
def plus(v,w):return[add(x,y)for x,y in zip(v,w)]
def vmul(v,a):return[scale(x,a)for x in v]
def dot(v,w):
 z=f()
 for r,(x,y)in enumerate(zip(v,w)):z=add(z,scale(mul(x,y),comb(N,r)))
 return z
def operator(v,which):
 if which=='z':return[scale(x,F(N,2)-r)for r,x in enumerate(v)]
 out=[f()for _ in range(D)]
 if which=='minus':
  for r in range(N):out[r+1]=scale(v[r],r+1)
 elif which=='plus':
  for r in range(1,D):out[r-1]=scale(v[r],N-r+1)
 return out
def tensors(v):
 z=operator(v,'z');lo=operator(v,'minus');hi=operator(v,'plus')
 return[z,plus(operator(z,'z'),vmul(v,-F(65,4))),lo,hi,
   plus(operator(lo,'z'),operator(z,'minus')),
   plus(operator(hi,'z'),operator(z,'plus')),
   operator(lo,'minus'),operator(hi,'plus')]
def rankmod(A,p=1000000007):
 A=[[(int(x.numerator)*pow(int(x.denominator),p-2,p))%p if isinstance(x,F) else int(x)%p for x in row]for row in A];rank=0
 for c in range(len(A[0])):
  pivot=next((r for r in range(rank,len(A))if A[r][c]),None)
  if pivot is None:continue
  A[rank],A[pivot]=A[pivot],A[rank];a=pow(A[rank][c],p-2,p);A[rank]=[(v*a)%p for v in A[rank]]
  for r in range(rank+1,len(A)):
   b=A[r][c]
   if b:A[r]=[(v-b*w)%p for v,w in zip(A[r],A[rank])]
  rank+=1
  if rank==len(A):break
 return rank
def primitive(row):
 den=lcm(*(x.denominator for x in row));out=[int(x*den)for x in row];g=gcd(*out)
 return[x//g for x in out]if g else out

def verify(root=None):
 start=time.time();root=Path(root)if root else Path(__file__).resolve().parent
 raw=(root/'certificate.json').read_bytes();c=json.loads(raw)
 assert(c['n'],c['t'])==(13,2)
 # Basis b_r=sqrt(binomial(13,r))D_r. u is unnormalized here.
 u=[f()for _ in range(D)];u[0]=f(0,5);u[9]=f(1)
 H=[[f()for _ in range(D)]for _ in range(7)]
 H[0][2]=f(0,1);H[0][11]=f(-1)
 for a,r in enumerate([3,4,5,6,12,13],1):H[a][r]=f(0,1)if r<=6 else f(1)
 gram=[[dot(v,w)for w in H]for v in H]
 assert all(z==f()for a,row in enumerate(gram)for b,z in enumerate(row)if a!=b)
 diag=[int(gram[a][a][0])for a in range(7)];assert diag==c['kernelGram']==[936,3146,7865,14157,18876,13,1]
 assert all(gram[a][a][1]==0 and diag[a]>0 for a in range(7))
 images=[u]+tensors(u);assert dot(u,u)==f(990)
 assert all(dot(u,w)==f()for w in images[1:])
 assert all(dot(h,w)==f()for h in H for w in images)
 # Rational 2x2 representation of field multiplication proves image rank>=7.
 doubled=[]
 for r in range(D):
  doubled.append([z for w in images for z in (w[r][0],11*w[r][1])])
  doubled.append([z for w in images for z in (w[r][1],w[r][0])])
 assert rankmod(doubled)==14
 A=[]
 for ell in range(8):
  mat=[]
  for a in range(7):
   row=[]
   for b in range(7):
    v=dot(H[a],tensors(H[b])[ell]);assert v[1]==0;row.append(v[0])
   mat.append(row)
  A.append(mat)
 assert A==[[[F(x)for x in row]for row in op]for op in c['compressedOperators']]
 pairs=list(combinations(range(7),2));charges=[2,3,4,5,6,3,4]
 sectors={q:[a for a,(i,j)in enumerate(pairs)if(charges[i]+charges[j])%9==q]for q in range(9)};sectors={q:v for q,v in sectors.items()if v}
 variables=[(q,a,b)for q,v in sectors.items()for a in range(len(v))for b in range(a,len(v))]
 assert pairs==list(map(tuple,c['pairs']))and charges==c['charges']
 assert sectors=={int(q):v for q,v in c['sectors'].items()}
 assert variables==list(map(tuple,c['variables']))and len(variables)==48
 def wedge(a):
  i,j=pairs[a];return[(i,j,1),(j,i,-1)]
 M=[[F(0)for _ in variables]for _ in range(8*49)]
 for col,(q,a,b)in enumerate(variables):
  aa,bb=sectors[q][a],sectors[q][b]
  for ca,cb in [(aa,bb)]+([(bb,aa)]if a!=b else[]):
   for x,p,sa in wedge(ca):
    for y,q2,sb in wedge(cb):
     for ell,op in enumerate(A):M[ell*49+p*7+q2][col]+=sa*sb*op[y][x]
 rows=[];meta=[]
 for i,row in enumerate(M):
  if any(row):rows.append(primitive(row));meta.append([i//49,(i%49)//7,i%7])
 assert rows==c['allPhysicalRows']and meta==c['allRowMetadata']and len(rows)==70
 selected=[rows[i]for i in c['selectedIndices']];assert selected==c['constraintRows']and len(selected)==40
 Z=c['kernel'];assert len(Z)==48 and all(len(row)==8 for row in Z)
 assert rankmod(selected)==40 and rankmod(Z)==8
 assert all(sum(row[i]*Z[i][j]for i in range(48))==0 for row in rows for j in range(8))
 diagonal=[i for i,(_,a,b)in enumerate(variables)if a==b];assert diagonal==c['diagonalVariableIndices']and len(diagonal)==21
 weights=c['positiveCoefficients'];assert len(weights)==21 and all(isinstance(x,int)and x>0 for x in weights)
 assert all(sum(weights[a]*Z[i][j]for a,i in enumerate(diagonal))==0 for j in range(8))
 report={'passed':True,'n':13,'order':2,'excludedDimension':3,'conditionalOn':'The code contains a collective rotation of u=sqrt(5/18)D0+sqrt(13/18)D9.',
 'method':'Independent exact quadratic-field spin actions and signed exterior contractions; Python standard library, no optimization or symbolic package.',
 'certificateSha256':hashlib.sha256(raw).hexdigest(),'tensorImageDimension':7,'orthogonalKernelDimension':7,'cyclicTwirlOrder':9,'exteriorDimension':21,'realInvariantCoordinates':48,'physicalContractionRows':70,'constraintRank':40,'kernelRank':8,'rankPrime':1000000007,'exactRowKernelProductsChecked':560,'strictlyPositiveDiagonalCoefficients':21,'exactPositiveFunctionalKernelProductsChecked':8,
 'scope':'All complex subspaces containing the specified state are excluded. This does not exclude arbitrary dimension-three subspaces.',
 'corollary':'No dimension-three second-order AC subspace contains a state supported on at most two Dicke components along any axis.',
 'elapsedSeconds':round(time.time()-start,3)}
 (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 import sys
 print(json.dumps(verify(sys.argv[1]if len(sys.argv)>1 else None),indent=2))
