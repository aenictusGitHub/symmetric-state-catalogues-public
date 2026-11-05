#!/usr/bin/env python3
"""Independent full ordered-word verification for a pure N15 AC4 exclusion.
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

def pullback(source,pptbasis,dual):
 start=time.time();n=source['n'];r=source['copies'];assert(n,r)==(15,4)
 by=defaultdict(list)
 for w in combinations_with_replacement(range(n+1),r):by[2*sum(w)-n*r].append(w)
 ix={M:{w:i for i,w in enumerate(ws)}for M,ws in by.items()}
 out=defaultdict(list)
 for i in range(n+1):
  for a in combinations_with_replacement(range(n+1),3):out[2*(i+sum(a))-4*n].append((i,a))
 expected={J:len(out[J])-len(out.get(J+2,[]))for J in range(0,4*n+1,2)if len(out[J])-len(out.get(J+2,[]))}
 assert {b['J2']:b['multiplicity']for b in pptbasis['blocks']}==expected
 assert sum((J+1)*p for J,p in expected.items())==(n+1)*comb(n+3,3)==13056
 Q={b['J2']:b['matrix']for b in dual['blocks']};precision=dual['precision'];power=precision+max(2*(max(map(int,b['dyadicColumnScales'])).bit_length()-1)for b in pptbasis['blocks'])
 W={M:[[0]*len(ws)for _ in ws]for M,ws in by.items()}
 raising=0;gram=0
 for block in pptbasis['blocks']:
  J=block['J2'];p=block['multiplicity'];H=[[int(x)for x in row]for row in block['highestWeight']];words=out[J]
  assert len(H)==len(words)and all(len(row)==p for row in H)
  # Reconstruct and test exact raising, in the averaged-triple basis.
  R=defaultdict(lambda:[0]*p)
  for (i,a),row in zip(words,H):
   if i<n:
    for k,x in enumerate(row):R[i+1,a][k]+=(i+1)*x
   for j,c in Counter(a).items():
    if j==n:continue
    aa=list(a);aa.remove(j);aa.append(j+1);target=(i,tuple(sorted(aa)))
    for k,x in enumerate(row):R[target][k]+=c*(j+1)*x
  assert all(x==0 for row in R.values()for x in row);raising+=p
  metric=[comb(n,i)*prod(comb(n,k)for k in a)*prod(factorial(c)for c in Counter(a).values())for i,a in words]
  for a in range(p):
   norm=sum(row[a]**2*g for row,g in zip(H,metric));assert norm>0
   if 'normNumerators'in block:assert norm==int(block['normNumerators'][a])
   for b in range(a):assert sum(row[a]*row[b]*g for row,g in zip(H,metric))==0;gram+=1
  scales=list(map(int,block['dyadicColumnScales']));assert all(s>0 and s&(s-1)==0 for s in scales);scale=max(scales);e=scale.bit_length()-1
  G=[[x*g*(scale//scales[a])for a,x in enumerate(row)]for row,g in zip(H,metric)]
  q=[[int(x)for x in row]for row in Q[J]];assert len(q)==p and all(len(row)==p for row in q)
  qcols=list(zip(*q));GQ=[[sum(x*y for x,y in zip(row,col))for col in qcols]for row in G]
  nz=[[(a,x)for a,x in enumerate(row)if x]for row in G]
  factor=(J+1)*2**(power-precision-2*e)
  for aidx,(i,a)in enumerate(words):
   ga=GQ[aidx]
   for bidx in range(aidx,len(words)):
    j,b=words[bidx];v=sum(ga[k]*x for k,x in nz[bidx])
    if not v:continue
    oi=n-j;pi=n-i;o=tuple(sorted(a+(oi,)));pw=tuple(sorted(b+(pi,)));M=2*sum(o)-n*r
    assert M==2*sum(pw)-n*r
    value=v*((-1)**(i+j))*(a.count(oi)+1)*(b.count(pi)+1)*factor
    u,vv=ix[M][o],ix[M][pw];W[M][u][vv]+=value
    if bidx!=aidx:W[M][vv][u]+=value
  print('independent output block',J,'time',time.time()-start,flush=True)
 numerators=[];matrices=[]
 for block in source['blocks']:
  J=block['J2'];p=block['multiplicity'];den=int(block['denominatorScale']);U=[[int(x)for x in row]for row in block['highestWeight']];D=[[0]*p for _ in range(p)]
  for step in range(J+1):
   M=J-2*step;Ucols=[[(i,row[a])for i,row in enumerate(U)if row[a]]for a in range(p)]
   WU=[[sum(row[i]*x for i,x in col)for col in Ucols]for row in W[M]]
   scale=den//comb(J,step);assert scale*comb(J,step)==den
   for a in range(p):
    for b in range(a,p):D[a][b]+=scale*sum(x*WU[i][b]for i,x in Ucols[a])
   if step==J:break
   nxt=[[0]*p for _ in by[M-2]]
   for iw,word in enumerate(by[M]):
    for i,c in Counter(word).items():
     if i==0:continue
     ww=list(word);ww.remove(i);ww.append(i-1);target=nxt[ix[M-2][tuple(sorted(ww))]];factor=c*(n-i+1)
     for a,x in enumerate(U[iw]):target[a]+=factor*x
   assert all(x%(step+1)==0 for row in nxt for x in row)
   U=[[x//(step+1)for x in row]for row in nxt]
  for a in range(p):
   for b in range(a,p):D[b][a]=D[a][b];numerators.append(D[a][b]*(1 if a==b else 2))
  matrices.append(D);print('independent input block',J,'time',time.time()-start,flush=True)
 return dict(denominator=576*2**power,numerators=numerators,matrixNumerators=matrices,outputHighestWeightsChecked=raising,outputCrossGramChecks=gram,seconds=time.time()-start)


def verify(directory=None):
 start=time.time();root=Path(directory or Path(__file__).resolve().parent);path=root/'certificate.json';raw=path.read_bytes();final=json.loads(raw);assert final['schema']=='ac4-four-copy-ppt-separation-v1';cert={'n':final['n'],'copies':final['copies'],'t':final['t'],'blocks':final['inputBlocks'],'variables':final['variables'],'selectedConstraints':final['constraints']};assert final['split']==[1,3]
 n,r,t=cert['n'],cert['copies'],cert['t'];assert(n,r,t)==(15,4,4);d=n+1
 words=list(combinations_with_replacement(range(d),r));by=defaultdict(list)
 for w in words:by[2*sum(w)-r*n].append(w)
 expected={M:len(by[M])-len(by.get(M+2,[]))for M in range(r*n%2,r*n+1,2)if len(by[M])-len(by.get(M+2,[]))}
 actual={b['J2']:b['multiplicity']for b in cert['blocks']};assert len(actual)==len(cert['blocks'])and actual==expected
 assert sum((M+1)*p for M,p in expected.items())==comb(n+r,r)==3876
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
 v=n*(n+2)
 polynomials=[[1 for i in range(d)],[(2*i-n)for i in range(d)],[3*(2*i-n)**2-v for i in range(d)],[(2*i-n)*(5*(2*i-n)**2-3*v+4)for i in range(d)],[35*(2*i-n)**4-(30*v-100)*(2*i-n)**2+3*v*(v-8)for i in range(d)]]
 assert all(sum(a*b for a,b in zip(polynomials[l],polynomials[h]))==0 for l in range(5)for h in range(l))
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
 # Independent rational pullback of the exact PPT functional.
 functional=pullback(cert,final['outputBasis'],final['outputDual'])
 f=[F(x,functional['denominator'])for x in functional['numerators']]
 weights_record=final['constraintWeights'];assert len(weights_record['numerators'])==len(A)==len(weights_record['denominatorExponents'])
 weights=[F(int(q),1<<int(e))if int(e)>=0 else F(int(q)*(1<<(-int(e))))for q,e in zip(weights_record['numerators'],weights_record['denominatorExponents'])]
 den=lcm(*(w.denominator for w in weights));wi=[w.numerator*(den//w.denominator)for w in weights]
 rowpart=[F(sum(w*row[i]for w,row in zip(wi,A)),den)for i in range(len(variables))]
 trace=[];scales=list(map(int,final['inputVariableScales']));assert len(scales)==len(variables)
 for iv,(ib,a,b)in enumerate(variables):
  J,p,dd,U=blocks[ib];ss=list(map(int,final['inputDyadicColumnScales'][ib]));assert len(ss)==p and all(s>0 and not s&(s-1)for s in ss)
  assert scales[iv]==dd*(J+1)*ss[a]*ss[b]
  trace.append(F(grams[ib][a][b]*dd*(J+1)*(1 if a==b else 2),factorial(r)**2))
 flat=[-tr-ff-rr for tr,ff,rr in zip(trace,f,rowpart)]
 Z=[[[F(0)for b in range(p)]for a in range(p)]for _,p,_,_ in blocks]
 for value,scale,(ib,a,b)in zip(flat,scales,variables):Z[ib][a][b]=Z[ib][b][a]=value/scale/(1 if a==b else 2)
 assert Z==[[[F(x)for x in row]for row in B]for B in final['basePositiveBlocks']]
 checks=[]
 for (J,p,_,_),B in zip(blocks,Z):checks.append({'J2':J,'dimension':p,'strictlyPositiveLDLPivots':[str(x)for x in rational_pd(B)]})
 pptchecks=[];expected_output={b['J2']:b['multiplicity']for b in final['outputBasis']['blocks']}
 assert len(expected_output)==len(final['outputDual']['blocks'])
 for B in final['outputDual']['blocks']:
  q=[[int(x)for x in row]for row in B['matrix']];assert len(q)==expected_output[B['J2']]
  pptchecks.append({'J2':B['J2'],'dimension':len(q),'strictlyPositiveLDLPivots':[str(x)for x in rational_pd(q)]})
 lowerw={0:F(9,88),4:F(35,88),11:F(35,88),15:F(9,88)}
 assert all(sum(w*(F(i)-F(15,2))**p for i,w in lowerw.items())==sum((F(i)-F(15,2))**p for i in range(16))/16 for p in range(4))
 report={'passed':True,'n':n,'t':t,'copies':r,'split':[1,3],'inputMomentDimension':len(words),'outputMomentDimension':13056,'variables':len(variables),'inputHighestWeightVectorsChecked':highestchecks,'physicalCrossGramChecks':normchecks,'selectedRowsReconstructed':len(A),'pptFunctionalCoordinatesReconstructed':len(f),'outputHighestWeightVectorsChecked':functional['outputHighestWeightsChecked'],'outputCrossGramChecks':functional['outputCrossGramChecks'],'certificateSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent full ordered-word reconstruction of the input basis and every selected first-leg contraction; exact 1|3 time-reversal pullback; rational separating identity; strict positive-definite LDL factorizations and their reconstructions. Standard library only, no discovery imports.','arbitraryComplexStatesCovered':True,'inputPositiveCertificate':{'blocks':len(checks),'positivePivots':sum(x['dimension']for x in checks),'allLDLReconstructionsExact':True,'certificates':checks},'pptPositiveCertificate':{'blocks':len(pptchecks),'positivePivots':sum(x['dimension']for x in pptchecks),'allLDLReconstructionsExact':True,'certificates':pptchecks},'exactDimensionOneThirdOrderWitness':{'weights':{str(i):str(w)for i,w in lowerw.items()},'momentsThroughOrder3Verified':True},'scope':'Universal exclusion of pure fourth-order AC states on 15 qubits; exact third-order witness independently checked.','statement':'K(15,t)=0 for every t>=4; maximum pure-state AC order is 3.','elapsedSeconds':time.time()-start}
 return report
if __name__=='__main__':
 root=Path(sys.argv[1])if len(sys.argv)>1 else Path(__file__).resolve().parent
 report=verify(root);(root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
