#!/usr/bin/env python3
"""Independent full ordered-word verification for a pure N25 AC6 source-identity component.
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
 start=time.time();root=Path(directory or Path(__file__).resolve().parent);path=root/'source-identity.json';raw=path.read_bytes();final=json.loads(raw);assert final['schema']=='exact-source-and-PPT-functional-identity-v1'
 cert={'n':final['n'],'copies':final['copies'],'t':final['t'],'blocks':final['inputBlocks'],'variables':final['variables'],'selectedConstraints':final['constraints']}
 n,r,t=cert['n'],cert['copies'],cert['t'];assert(n,r,t)==(25,4,6);d=n+1
 words=list(combinations_with_replacement(range(d),r));by=defaultdict(list)
 for w in words:by[2*sum(w)-r*n].append(w)
 expected={M:len(by[M])-len(by.get(M+2,[]))for M in range(r*n%2,r*n+1,2)if len(by[M])-len(by.get(M+2,[]))}
 actual={b['J2']:b['multiplicity']for b in cert['blocks']};assert len(actual)==len(cert['blocks'])and actual==expected
 assert sum((M+1)*p for M,p in expected.items())==comb(n+r,r)==23751
 # Store symmetric ordered amplitudes by permutation orbit. Orbit sizes in
 # the metric reproduce the full ordered tensor-product inner product exactly.
 metric={w:factorial(r)//prod(factorial(c)for c in Counter(w).values())*prod(comb(n,i)for i in w)for w in words}
 prime=1000003;assert all(prime%q for q in range(2,isqrt(prime)+1))
 def step(vector,direction,divisor=1):
  out=defaultdict(int)
  for w,v in vector.items():
   counts=Counter(w)
   for i in counts:
    nxt=i+direction
    if not 0<=nxt<=n:continue
    target=list(w);target.remove(i);target.append(nxt);target=tuple(sorted(target))
    # Coefficient of an ordered target receives one preimage per target slot.
    out[target]+=(counts.get(nxt,0)+1)*(i+1 if direction==1 else n-i+1)*v
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
    assert word not in v;v[word]=factor*row[a]
   assert not step(v,1);vectors.append(v);highestchecks+=1
  G=[[inner(a,b)for b in vectors]for a in vectors];assert all(G[a][a]>0 for a in range(p));assert all(G[a][b]==0 for a in range(p)for b in range(a));grams.append(G);U={J:vectors}
  for k in range(1,J+1):
   vectors=[step(v,-1,k)for v in vectors];U[J-2*k]=vectors
   for a in range(p):
    for b in range(a,p):assert inner(vectors[a],vectors[b])==comb(J,k)*G[a][b];normchecks+=1
  den=int(record['denominatorScale']);assert den>0 and all(den%comb(J,k)==0 for k in range(J+1))
  blocks.append((J,p,den,U));print('physical block',ib,'seconds',time.time()-start,flush=True)
 variables=[[ib,a,b]for ib,(_,p,_,_)in enumerate(blocks)for a in range(p)for b in range(a,p)];assert variables==cert['variables']
 # Positive physical control: the coherent-state twirl is the normalized
 # highest-spin projector. Its ordered coefficients have a closed Vandermonde
 # form, independent of the certificate's basis or contraction row selection.
 topib=max(range(len(blocks)),key=lambda k:blocks[k][0]);topJ,topp,topden,topU=blocks[topib]
 assert topJ==r*n and topp==1
 topamplitude=topU[topJ][0][(n,)*r]
 for M,vec in topU.items():
  assert set(vec[0])==set(by[M]) and set(vec[0].values())=={topamplitude}
  assert sum(metric[w]for w in by[M])==comb(r*n,(r*n-M)//2)
 coherentB=F(factorial(r)**2,topden*(topJ+1)*topamplitude**2)
 topiv=variables.index([topib,0,0]);coherentchecks=0;coherentnonzero=0

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
  for J,p,den,U in blocks:
   pairs=[(a,b)for a in range(p)for b in range(a,p)];values=[0]*len(pairs)
   for i in range(d):
    m=M+2*i-n
    if m not in U:continue
    leftword=tuple(sorted((i,)+left));rightword=tuple(sorted((i,)+right))
    lv=[v.get(leftword,0)for v in U[m]];rv=[v.get(rightword,0)for v in U[m]]
    if not any(lv)or not any(rv):continue
    scale=(den//comb(J,(J-m)//2))*comb(n,i)*polynomials[ell][i]
    if not scale:continue
    for k,(a,b)in enumerate(pairs):values[k]+=scale*(lv[a]*rv[b]+(lv[b]*rv[a]if a!=b else 0))
   row.extend(values)
  closed=sum(F(comb(n,i)*polynomials[ell][i],(r*n+1)*comb(r*n,(r*n-(M+2*i-n))//2))for i in range(d))
  assert F(row[topiv])*coherentB/factorial(r)**2==closed
  coherentchecks+=1;coherentnonzero+=int(closed!=0)
  common=gcd(*row);assert common;row=[x//common for x in row]
  if 'selectedRows'in cert:assert row==[int(x)for x in cert['selectedRows'][ri]]
  A.append(row)
  if ri%100==0:print('physical row',ri,'seconds',time.time()-start,flush=True)
 # Reconstruct the physical trace from every magnetic descendant, with (4!)^2 normalization.
 trace=[F(0)]*len(variables);traceschecked=0
 for ib,(J,p,dd,U)in enumerate(blocks):
  G=[[F(0)for _ in range(p)]for _ in range(p)]
  for M,vectors in U.items():
   k=(J-M)//2
   for a in range(p):
    for b in range(a,p):G[a][b]+=F(dd*inner(vectors[a],vectors[b]),comb(J,k)*factorial(r)**2)
  for a in range(p):
   for b in range(a,p):assert G[a][b]==F(dd*(J+1)*grams[ib][a][b],factorial(r)**2)
  for iv,(bi,a,b)in enumerate(variables):
   if bi==ib:trace[iv]=G[a][b]*(1 if a==b else 2)
  traceschecked+=1
 weights_record=final['constraintWeights'];assert len(weights_record['numerators'])==len(A)==len(weights_record['denominatorExponents'])
 weights=[F(int(q),1<<int(e))if int(e)>=0 else F(int(q)*(1<<(-int(e))))for q,e in zip(weights_record['numerators'],weights_record['denominatorExponents'])]
 den=lcm(*(w.denominator for w in weights));wi=[w.numerator*(den//w.denominator)for w in weights]
 rowpart=[F(sum(w*row[i]for w,row in zip(wi,A)),den)for i in range(len(variables))]
 scales=list(map(int,final['inputVariableScales']));assert len(scales)==len(variables)
 for iv,(ib,a,b)in enumerate(variables):
  J,p,dd,U=blocks[ib];ss=list(map(int,final['inputDyadicColumnScales'][ib]));assert len(ss)==p and all(s>0 and not s&(s-1)for s in ss)
  assert scales[iv]==dd*(J+1)*ss[a]*ss[b]
 ppt=final['pptFunctional'];pf=[F(int(x),int(ppt['denominator']))for x in ppt['numerators']];assert len(pf)==len(trace)
 flat=[-tr-rr-ff for tr,rr,ff in zip(trace,rowpart,pf)]
 Z=[[[F(0)for b in range(p)]for a in range(p)]for _,p,_,_ in blocks]
 for value,scale,(ib,a,b)in zip(flat,scales,variables):Z[ib][a][b]=Z[ib][b][a]=value/scale/(1 if a==b else 2)
 assert Z==[[[F(x)for x in row]for row in B]for B in final['basePositiveBlocks']]
 checks=[]
 for (J,p,_,_),B in zip(blocks,Z):checks.append({'J2':J,'dimension':p,'strictlyPositiveLDLPivots':[str(x)for x in rational_pd(B)]})
 report={'passed':True,'n':n,'t':t,'copies':r,'inputMomentDimension':len(words),'variables':len(variables),'inputHighestWeightVectorsChecked':highestchecks,'physicalCrossGramChecks':normchecks,'selectedRowsReconstructed':len(A),'allDescendantTraceBlocksReconstructed':traceschecked,'coherentStateControl':{'passed':True,'closedVandermondeRowsChecked':coherentchecks,'nonzeroRowsConfirmNonAnticoherence':coherentnonzero,'allOrderedDescendantAmplitudesConstant':True,'normalizedHighestSpinProjectorTrace':'1'},'certificateSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent ordered-word reconstruction, stored by exact permutation orbits, of every input block and selected first-leg contraction equation; trace independently summed over all magnetic descendants; exact rational separating identity; strictly positive LDL factorizations with reconstruction. Standard library only, no discovery imports.','arbitraryComplexStatesCovered':True,'inputPositiveCertificate':{'blocks':len(checks),'positivePivots':sum(x['dimension']for x in checks),'allLDLReconstructionsExact':True,'certificates':checks},'sourceIdentityOnly':True,'pptFunctionalIndependentlyVerified':False,'scope':'Exact physical source-identity component for N25 AC6, treating the saved rational PPT functional as an algebraic input. Nonnegativity of that functional needs a separate physical PPT verification.','statement':'Source trace and contraction rows reconstructed; exact source identity and source positivity verified. This report alone is not a global exclusion proof.','elapsedSeconds':time.time()-start}
 return report
if __name__=='__main__':
 root=Path(sys.argv[1])if len(sys.argv)>1 else Path(__file__).resolve().parent
 report=verify(root);(root/'source-identity-independent.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
