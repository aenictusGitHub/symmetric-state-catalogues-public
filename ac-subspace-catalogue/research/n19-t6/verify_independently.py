#!/usr/bin/env python3
"""Independent full ordered-word verification for a pure N19 AC6 exclusion.
Only Python standard library arithmetic. No discovery code is imported.
"""
from pathlib import Path
from itertools import combinations_with_replacement,permutations,product
from collections import Counter,defaultdict
from math import comb,prod,factorial,gcd,lcm,isqrt
from fractions import Fraction as F
import json,sys,time,hashlib

def modular_rank(matrix,prime):
 a=[[int(x)%prime for x in row]for row in matrix];r=0
 for col in range(len(a[0])):
  pivot=next((i for i in range(r,len(a))if a[i][col]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];inv=pow(a[r][col],-1,prime);a[r]=[(x*inv)%prime for x in a[r]]
  for i in range(r+1,len(a)):
   if a[i][col]:
    c=a[i][col];a[i]=[(x-c*y)%prime for x,y in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r

def rational_pd(matrix):
 a=[[F(x)for x in row]for row in matrix];n=len(a);L=[[F(i==j)for j in range(n)]for i in range(n)];D=[]
 assert n and all(len(row)==n for row in a) and all(a[i][j]==a[j][i]for i in range(n)for j in range(n))
 for j in range(n):
  p=a[j][j]-sum(L[j][h]**2*D[h]for h in range(j));assert p>0,(j,p);D.append(p)
  for i in range(j+1,n):L[i][j]=(a[i][j]-sum(L[i][h]*L[j][h]*D[h]for h in range(j)))/p
 for i in range(n):
  for j in range(i+1):assert sum(L[i][h]*D[h]*L[j][h]for h in range(min(i,j)+1))==a[i][j]
 return D


def verify(directory=None):
 start=time.time();root=Path(directory or Path(__file__).resolve().parent);path=root/'certificate.json';raw=path.read_bytes();cert=json.loads(raw);assert cert['schema']=='three-copy-ac6-positive-probability-v1'
 n,r,t=cert['n'],cert['copies'],cert['t'];assert(n,r,t)==(19,3,6);d=n+1
 words=list(combinations_with_replacement(range(d),r));by=defaultdict(list)
 for w in words:by[2*sum(w)-r*n].append(w)
 expected={M:len(by[M])-len(by.get(M+2,[]))for M in range(r*n%2,r*n+1,2)if len(by[M])-len(by.get(M+2,[]))}
 actual={b['J2']:b['multiplicity']for b in cert['blocks']};assert len(actual)==len(cert['blocks'])and actual==expected
 assert sum((M+1)*p for M,p in expected.items())==comb(n+r,r)==1540
 metric={w:prod(comb(n,i)for i in w)for w in product(range(d),repeat=r)}
 prime=1000003;assert all(prime%q for q in range(2,isqrt(prime)+1))
 def step(vector,direction,divisor=1):
  out=defaultdict(int)
  for w,v in vector.items():
   for k,i in enumerate(w):
    nxt=i+direction
    if not 0<=nxt<=n:continue
    target=list(w);target[k]=nxt
    out[tuple(target)]+=(i+1 if direction==1 else n-i+1)*v
  result={}
  for w,v in out.items():
   assert v%divisor==0
   if v:result[w]=v//divisor
  return result
 def inner(a,b):return sum(v*b.get(w,0)*metric[w]for w,v in a.items())
 blocks=[];grams=[];normchecks=0;highestchecks=0
 for ib,record in enumerate(cert['blocks']):
  J=record['J2'];p=record['multiplicity'];H=[[int(x)for x in row]for row in record['highestWeight']]
  assert len(H)==len(by[J])and all(len(row)==p for row in H);assert modular_rank(H,prime)==p
  vectors=[]
  for a in range(p):
   v={}
   for word,row in zip(by[J],H):
    if not row[a]:continue
    # r! times an averaged occupation word, independently expanded.
    factor=prod(factorial(c)for c in Counter(word).values())
    for ordered in set(permutations(word)):
     assert ordered not in v;v[ordered]=factor*row[a]
   assert not step(v,1);vectors.append(v);highestchecks+=1
  G=[[inner(a,b)for b in vectors]for a in vectors];assert all(G[a][a]>0 for a in range(p));assert all(G[a][b]==0 for a in range(p)for b in range(a));grams.append(G);U={J:vectors}
  for k in range(1,J+1):
   vectors=[step(v,-1,k)for v in vectors];U[J-2*k]=vectors
   for a in range(p):
    for b in range(a,p):assert inner(vectors[a],vectors[b])==comb(J,k)*G[a][b];normchecks+=1
  den=int(record['denominatorScale']);assert den>0 and all(den%comb(J,k)==0 for k in range(J+1))
  blocks.append((J,p,den,U));print('physical block',ib,'seconds',time.time()-start,flush=True)
 variables=[[ib,a,b]for ib,(_,p,_,_)in enumerate(blocks)for a in range(p)for b in range(a,p)];assert variables==cert['variables']
 lower=defaultdict(list)
 for w in combinations_with_replacement(range(d),r-1):lower[2*sum(w)-(r-1)*n].append(w)
 polynomials=[]
 for ell in range(t+1):
  v=[F((2*i-n)**ell)for i in range(d)]
  for u in polynomials:
   factor=sum(a*b for a,b in zip(v,u))/sum(a*a for a in u)
   v=[a-factor*b for a,b in zip(v,u)]
  den=lcm(*(x.denominator for x in v));v=[int(x*den)for x in v];common=gcd(*v)
  polynomials.append([x//common for x in v])
 assert all(sum(a*b for a,b in zip(polynomials[l],polynomials[h]))==0 for l in range(t+1)for h in range(l))
 A=[]
 for ri,(ell,M,ia,ib)in enumerate(cert['selectedConstraints']):
  assert 1<=ell<=t and 0<=ia<=ib<len(lower[M]);left,right=lower[M][ia],lower[M][ib];row=[]
  for blockid,a,b in variables:
   J,_,den,U=blocks[blockid];value=0
   for i in range(d):
    m=M+2*i-n
    if m not in U:continue
    p=U[m][a].get((i,)+left,0)*U[m][b].get((i,)+right,0)
    if a!=b:p+=U[m][b].get((i,)+left,0)*U[m][a].get((i,)+right,0)
    value+=(den//comb(J,(J-m)//2))*comb(n,i)*polynomials[ell][i]*p
   row.append(value)
  common=gcd(*row);assert common;row=[x//common for x in row]
  if 'selectedRows'in cert:assert row==[int(x)for x in cert['selectedRows'][ri]]
  A.append(row)
  if ri%100==0:print('physical row',ri,'seconds',time.time()-start,flush=True)
 # Verify a complete exact kernel of all selected necessary equations.
 rank=modular_rank(A,prime);assert rank==len(A)==cert['rank']==88
 kernel=[[int(x)for x in row]for row in cert['kernel']];dimension=len(variables)-rank
 assert dimension==6 and len(kernel)==len(variables) and all(len(row)==dimension for row in kernel)
 assert modular_rank(kernel,prime)==dimension
 assert all(sum(a*b for a,b in zip(row,col))==0 for row in A for col in zip(*kernel))
 # Reconstruct every physical total-spin probability by summing all descendants.
 probabilityrows=[];wholetrace=[F(0)]*len(variables);totalspinchecks=0
 for ib,(J,p,dd,U)in enumerate(blocks):
  blocktrace=[[F(0)for b in range(p)]for a in range(p)]
  for M,vectors in U.items():
   k=(J-M)//2
   for a in range(p):
    for b in range(a,p):blocktrace[a][b]+=F(dd*inner(vectors[a],vectors[b]),comb(J,k)*factorial(r)**2)
  expected=[[F(dd*(J+1)*grams[ib][a][b],factorial(r)**2)for b in range(p)]for a in range(p)]
  assert all(blocktrace[a][b]==expected[a][b]for a in range(p)for b in range(a,p))
  row=[F(0)]*len(variables)
  for iv,(bi,a,b)in enumerate(variables):
   if bi==ib:row[iv]=blocktrace[a][b]*(1 if a==b else 2)
  wholetrace=[a+b for a,b in zip(wholetrace,row)];probabilityrows.append(row);totalspinchecks+=1
 tracekernel=[sum(a*b for a,b in zip(wholetrace,col))for col in zip(*kernel)]
 probability=cert['probabilityCertificate'];indices=probability['blockIndices'];spins=probability['spins2']
 assert len(indices)==len(spins)==6 and len(set(indices))==6
 assert [blocks[i][0]for i in indices]==spins
 forms=[[sum(a*b for a,b in zip(probabilityrows[i],col))for col in zip(*kernel)]for i in indices]
 assert forms==[[F(x)for x in row]for row in probability['kernelRows']]
 assert tracekernel==[F(x)for x in probability['traceKernelRow']]
 weight=F(probability['traceWeight']);weights=[F(x)for x in probability['positiveWeights']]
 assert weight>0 and all(x>0 for x in weights)
 assert all(weight*tr+sum(w*row[c]for w,row in zip(weights,forms))==0 for c,tr in enumerate(tracekernel))
 report={'passed':True,'n':n,'t':t,'copies':r,'inputMomentDimension':len(words),'realInvariantVariables':len(variables),'constraintRank':rank,'kernelDimension':dimension,
 'inputHighestWeightVectorsChecked':highestchecks,'physicalCrossGramChecks':normchecks,'selectedRowsReconstructed':len(A),'totalSpinProbabilitiesReconstructed':totalspinchecks,
 'certificateSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent full ordered-word reconstruction of the complete input representation, every selected contraction row, exact complete kernel and all physical total-spin traces; rational positive trace-versus-probabilities identity. Standard library only; no discovery imports.',
 'arbitraryComplexStatesCovered':True,
 'probabilityCertificate':{'spins2':spins,'traceWeight':str(weight),'strictlyPositiveWeights':[str(w)for w in weights],'allDescendantTraceReconstructionsExact':True,'kernelIdentityExact':True,'contradiction':'Trace one and nonnegative total-spin probabilities make the certified expression strictly positive, but the necessary AC6 equations force it to zero.'},
 'scope':'Universal exclusion of pure sixth-order AC states on 19 qubits.','statement':'K(19,t)=0 for every t>=6.','elapsedSeconds':time.time()-start}
 return report
if __name__=='__main__':
 root=Path(sys.argv[1])if len(sys.argv)>1 else Path(__file__).resolve().parent
 report=verify(root);(root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
