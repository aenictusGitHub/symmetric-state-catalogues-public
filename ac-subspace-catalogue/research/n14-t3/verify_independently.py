#!/usr/bin/env python3
"""Independent physical-tensor-word verification of the N14 AC3 Young22 obstruction.
No discovery code is imported. Uses only the Python 3 standard library, with exact
integer and rational arithmetic; all physical reconstruction is separate.
"""
from pathlib import Path
from itertools import combinations,product
from collections import defaultdict,Counter
from math import comb,gcd,prod,isqrt,lcm
import json,time,hashlib,sys
from fractions import Fraction as F

def modular_rank(matrix,prime):
 a=[[int(x)%prime for x in row]for row in matrix];r=0
 for col in range(len(a[0])):
  pivot=next((i for i in range(r,len(a))if a[i][col]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];inverse=pow(a[r][col],-1,prime)
  a[r]=[(x*inverse)%prime for x in a[r]]
  for i in range(r+1,len(a)):
   if a[i][col]:
    scale=a[i][col];a[i]=[(x-scale*y)%prime for x,y in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r

def rational_pd(matrix):
 # Fraction LDL without square roots: positive pivots certify positive definiteness.
 a=[[F(x) for x in row] for row in matrix];n=len(a)
 assert n>0 and all(len(row)==n for row in a)
 assert all(a[i][j]==a[j][i] for i in range(n)for j in range(n))
 lower=[[F(int(i==j))for j in range(n)]for i in range(n)];pivots=[]
 for j in range(n):
  pivot=a[j][j]-sum(lower[j][h]**2*pivots[h]for h in range(j))
  assert pivot>0, f'Nonpositive LDL pivot {j}'
  pivots.append(pivot)
  for i in range(j+1,n):
   lower[i][j]=(a[i][j]-sum(lower[i][h]*lower[j][h]*pivots[h]for h in range(j)))/pivot
 # Reconstruct the original symmetric matrix to audit every entry.
 for i in range(n):
  for j in range(i+1):assert sum(lower[i][h]*pivots[h]*lower[j][h]for h in range(min(i,j)+1))==a[i][j]
 return pivots

def audit_row_dual(A,multipliers,blockrecords,variables,savedblocks=None):
 weights=[F(x)for x in multipliers];assert len(weights)==len(A)
 # Exact combination of the independently reconstructed necessary equations.
 denominator=lcm(*(x.denominator for x in weights))
 integer_weights=[x.numerator*(denominator//x.denominator)for x in weights]
 flat=[F(sum(w*row[i]for w,row in zip(integer_weights,A)),denominator)for i in range(len(variables))]
 blocks=[[[F(0)for _ in range(b['multiplicity'])]for _ in range(b['multiplicity'])]for b in blockrecords]
 for value,(ib,a,b)in zip(flat,variables):
  blocks[ib][a][b]=value if a==b else value/2
  blocks[ib][b][a]=blocks[ib][a][b]
 if savedblocks is not None:
  assert blocks==[[[F(x)for x in row]for row in b]for b in savedblocks]
 certificates=[]
 for b,record in zip(blocks,blockrecords):
  pivots=rational_pd(b)
  certificates.append({'J2':record['J2'],'dimension':len(b),'strictlyPositiveLDLPivots':[str(x)for x in pivots]})
 return {'rowCombinationOfNecessaryEquationsReconstructedExactly':True,
         'multiplierCount':len(weights),'strictlyPositiveDefiniteDualBlocks':len(blocks),
         'allLDLReconstructionsExact':True,'positiveDefinitenessCertificates':certificates,
         'conclusion':'The exact linear combination of necessary equations is a strictly positive functional on every nonzero positive semidefinite invariant moment. Thus the normalized moment cannot exist.'}

def verify(directory=None):
 root=Path(directory) if directory is not None else Path(__file__).resolve().parent
 path=root/'certificate.json';raw=path.read_bytes();data=json.loads(raw);n=data['n'];d=n+1;start=time.time()
 assert(n,data['t'],data['excludedDimension'])==(14,3,2)
 assert data['schurPartition']==[2,2] and data['momentDimension']==4200
 polynomials=[[1 for i in range(d)],[(2*i-n)for i in range(d)],
              [3*(2*i-n)**2-n*(n+2)for i in range(d)],
              [(2*i-n)*(5*(2*i-n)**2-3*n*(n+2)+4)for i in range(d)]]
 assert all(sum(x*y for x,y in zip(polynomials[a],polynomials[b]))==0 for a in range(4)for b in range(a))
 prime=int(data.get('determinantModulus',1000003));assert prime>2 and all(prime%p for p in range(2,isqrt(prime)+1))
 def rank(A):return modular_rank(A,prime)
 pairs=list(combinations(range(d),2));words=[p+q for p in pairs for q in pairs if p<=q]
 by=defaultdict(list)
 for w in words:by[2*sum(w)-4*n].append(w)
 weights={w:prod(comb(n,r)for r in w)for w in product(range(d),repeat=4)}
 def orbit(w):
  a,b,c,e=w;out={}
  for left,right in [((a,b),(c,e)),((c,e),(a,b))]:
   for p,sp in [(left,1),(left[::-1],-1)]:
    for q,sq in [(right,1),(right[::-1],-1)]:
     word=p+q;sign=sp*sq
     assert word not in out or out[word]==sign
     out[word]=sign
  return out
 def step(vector,direction,divisor=1):
  out=defaultdict(int)
  for source,value in vector.items():
   for leg,r in enumerate(source):
    target=list(source);target[leg]+=direction
    if not 0<=target[leg]<=n:continue
    factor=n-r+1 if direction==-1 else r+1
    out[tuple(target)]+=factor*value
  result={}
  for word,value in out.items():
   assert value%divisor==0
   if value:result[word]=value//divisor
  return result
 def inner(u,v):return sum(value*v.get(w,0)*weights[w]for w,value in u.items())
 quad=Counter(2*sum(q)-4*n for q in combinations(range(d),4))
 weightdim={M:len(by.get(M,[]))-quad[M]for M in range(0,4*n+3,2)}
 expected={M:weightdim[M]-weightdim.get(M+2,0)for M in range(0,4*n+1,2)if weightdim[M]-weightdim.get(M+2,0)}
 actual={b['J2']:b['multiplicity']for b in data['blocks']}
 assert len(actual)==len(data['blocks'])and actual==expected
 assert sum((M+1)*m for M,m in expected.items())==d*d*(d*d-1)//12==4200
 blocks=[];highest_checks=0;norm_checks=0
 for ib,block in enumerate(data['blocks']):
  M=block['J2'];mult=block['multiplicity'];H=[[int(x)for x in row]for row in block['highestWeight']];physical=[]
  assert len(H)==len(by[M])and all(len(row)==mult for row in H)
  assert rank(H)==mult
  for a in range(mult):
   vec={}
   for word,row in zip(by[M],H):
    if row[a]:
     for w,sign in orbit(word).items():
      assert w not in vec;vec[w]=sign*row[a]
   assert not step(vec,1)
   for q in combinations(range(d),4):
    if 2*sum(q)-4*n==M:
     a0,b0,c0,e0=q
     assert vec.get((a0,b0,c0,e0),0)-vec.get((a0,c0,b0,e0),0)+vec.get((a0,e0,b0,c0),0)==0
   physical.append(vec);highest_checks+=1
  gram=[[inner(a,b)for b in physical]for a in physical]
  assert all(gram[a][b]==0 for a in range(mult)for b in range(a))
  U={M:physical}
  for k in range(1,M+1):
   physical=[step(v,-1,k)for v in physical];U[M-2*k]=physical
   for a in range(mult):
    for b in range(a,mult):
     assert inner(physical[a],physical[b])==comb(M,k)*gram[a][b];norm_checks+=1
  den=int(block['denominatorScale']);assert den>0 and all(den%comb(M,k)==0 for k in range(M+1))
  blocks.append({'J2':M,'U':U,'scale':den})
  print('physical block',ib,'highest',highest_checks,'seconds',time.time()-start,flush=True)
 variables=[[ib,a,b]for ib,B in enumerate(data['blocks'])for a in range(B['multiplicity'])for b in range(a,B['multiplicity'])]
 assert data['variables']==variables and len(variables)==1055
 low=defaultdict(list)
 for a in range(d):
  for b,c in pairs:low[2*(a+b+c)-3*n].append((a,b,c))
 A=[]
 for rowindex,(ell,M3,alpha,beta)in enumerate(data['selectedConstraints']):
  assert 1<=ell<=3 and 0<=alpha<=beta<len(low[M3])
  left=low[M3][alpha];right=low[M3][beta];row=[]
  for ib,a,b in variables:
   B=blocks[ib];J=B['J2'];value=0
   for i in range(d):
    M=M3+2*i-n
    if M not in B['U']:continue
    v=B['U'][M];z=2*i-n;f=z if ell==1 else 3*z*z-n*(n+2)if ell==2 else z*(5*z*z-3*n*(n+2)+4)
    xy=v[a].get((i,)+left,0)*v[b].get((i,)+right,0)
    if a!=b:xy+=v[b].get((i,)+left,0)*v[a].get((i,)+right,0)
    value+=(B['scale']//comb(J,(J-M)//2))*comb(n,i)*f*xy
   row.append(value)
  common=gcd(*row);row=[x//common for x in row] if common else row
  if 'selectedRows' in data:assert row==[int(x)for x in data['selectedRows'][rowindex]]
  A.append(row)
  if rowindex%100==0:print('physical rows',rowindex+1,'seconds',time.time()-start,flush=True)
 assert len(A)==len(data['selectedConstraints'])
 report={'passed':True,'n':n,'t':3,'excludedDimension':2,'certificateSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent full physical-word expansion, physical spin-ladder operations and direct first-leg contractions; no canonical ladder or contraction code imported. The necessary integer rows are reconstructed exactly; no rank computation of the constraint matrix is needed.','young22Dimension':4200,'schurMultiplicitiesIndependentlyCounted':True,'highestWeightVectorsIndependentlyChecked':highest_checks,'loweringCrossGramIdentitiesIndependentlyChecked':norm_checks,'constraintRowsReconstructed':len(A),'variables':len(variables),'highestWeightRankModulus':prime,'physicalConstraintRowsIndependentlyReconstructed':True,'noFloatingPointProofSteps':True,'scope':'Independent reconstruction and modular rank audit only; this report alone is not a nonexistence proof.','elapsedSeconds':time.time()-start}
 assert 'rowMultipliers' in data
 if 'rowMultipliers' in data:
  report['positiveDefiniteRowCertificate']=audit_row_dual(A,data['rowMultipliers'],data['blocks'],variables)
  report['method']='Independent full physical-word reconstruction of all highest weights, spin ladders and necessary contraction rows; exact rational row combination and positive-definite LDL certificate.'
  report['scope']='Universal exclusion of dimension-two third-order anticoherent subspaces of spin7, including arbitrary complex bases.'
  lower={0:F(4,21),7:F(13,21),14:F(4,21)}
  assert all(sum(w*F(i-7)**p for i,w in lower.items())==sum(F(i-7)**p for i in range(15))/15 for p in range(4))
  assert all(abs(i-j)>3 for i in lower for j in lower if i!=j)
  report['exactDimensionOneWitness']={'support':[0,7,14],'squaredAmplitudes':['4/21','13/21','4/21'],'diagonalMomentsThroughRankThreeExact':True,'nonzeroMagneticComponentsVanishBySupportGaps':True}
  report['statement']='K(14,3)=1: no dimension-two third-order AC subspace exists, and an exact three-component third-order state supplies dimension one.'
  report['elapsedSeconds']=time.time()-start
 return report
if __name__=='__main__':
 root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
 report=verify(root);out=root/'independent-verification.json'
 out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
