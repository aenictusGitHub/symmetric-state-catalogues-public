#!/usr/bin/env python3
"""Physical exact PPT-functional verification; Python standard library only.
Reconstructs both rational coupled bases and every raw-source coefficient.
No numerical discovery data, matrix library, tolerance, or saved overlap is used.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from collections import Counter,defaultdict
from math import comb,factorial,prod,lcm,gcd
import json,time,hashlib,sys
from coupled_engine import cg_highest,coupled_vector,lower_vector,sixj_raw

def verify(directory=None):
 root=Path(directory or Path(__file__).resolve().parent);start=time.time()
 sourceRaw=(root/'source-identity.json').read_bytes();source=json.loads(sourceRaw)
 cutsRaw=(root/'rational-cuts.json').read_bytes();data=json.loads(cutsRaw)
 assert(source['n'],source['copies'],source['t'])==(25,4,6) and(data['n'],data['copies'])==(25,4)
 n=25;BY={}
 for r in [1,2,3,4]:
  by=defaultdict(list)
  for word in combinations_with_replacement(range(26),r):by[2*sum(word)-r*n].append(word)
  BY[r]=by
 def gnum(word):return prod(comb(n,i)**v*factorial(v) for i,v in Counter(word).items())
 reps={('one',25,0):[{(25-q,):1} for q in range(26)]};norms={('one',25,0):Q(1)}
 parentChecks=0;descChecks=0
 for name,r,records in [('pair',2,data['parentPairBlocks']),('triple',3,data['parentTripleBlocks'])]:
  expected={M:len(BY[r][M])-len(BY[r].get(M+2,[])) for M in range(r*n%2,r*n+1,2) if len(BY[r][M])-len(BY[r].get(M+2,[]))}
  assert {B['J2']:B['multiplicity'] for B in records}==expected
  assert sum((S+1)*p for S,p in expected.items())==comb(n+r,r)
  for B in records:
   S=B['J2'];H=[[int(v) for v in row] for row in B['highestWeight']];p=B['multiplicity'];ws=BY[r][S]
   assert len(H)==len(ws) and all(len(row)==p for row in H)
   for a in range(p):
    v={w:row[a] for w,row in zip(ws,H) if row[a]};raising={}
    for w,value in v.items():
     for i,count in Counter(w).items():
      if i==n:continue
      z=list(w);z.remove(i);z.append(i+1);z=tuple(sorted(z));raising[z]=raising.get(z,0)+count*(i+1)*value
    assert all(value==0 for value in raising.values())
    norm=sum(Q(value*value*gnum(w),factorial(r)) for w,value in v.items());assert norm>0
    for b in range(a):assert sum(Q(row[a]*row[b]*gnum(w),factorial(r)) for w,row in zip(ws,H))==0
    key=(name,S,a);norms[key]=norm;desc=[v]
    for q in range(1,S+1):desc.append(lower_vector(n,desc[-1],q))
    for q,u in enumerate(desc):assert sum(Q(value*value*gnum(w),factorial(r)) for w,value in u.items())==comb(S,q)*norm;descChecks+=1
    reps[key]=desc;parentChecks+=1
 def normcouple(a,b,J):return Q(cg_highest(a[1],b[1],J)[3])*norms[a]*norms[b]
 def symcouple(a,b,J):
  W,_=coupled_vector(a[1],reps[a],b[1],reps[b],J);out={}
  for (u,v),value in W.items():w=tuple(sorted(u+v));out[w]=out.get(w,0)+value
  return {w:value for w,value in out.items() if value}
 def channels(split,J):
  if split=='22':return [(('pair',A,0),('pair',B,0)) for A in range(2,51,4) for B in range(2,51,4) if abs(A-B)<=J<=A+B]
  return [(('one',25,0),('triple',B['J2'],a)) for B in data['parentTripleBlocks'] if abs(25-B['J2'])<=J<=25+B['J2'] for a in range(B['multiplicity'])]
 # Every cut is a positive rational weight times a highest-vector projector.
 QRAW={};LAY={};cutIds=set();splitCounts=defaultdict(int)
 for cut in data['cuts']:
  assert cut['cutId'] not in cutIds;cutIds.add(cut['cutId']);split='22' if cut['split']==[2,2] else '13';assert cut['split'] in [[1,3],[2,2]]
  L=cut['J2'];keys=[]
  for label in cut['layout']:
   if split=='22':keys.append((('pair',label[0],0),('pair',label[1],0)))
   else:keys.append((('one',25,0),('triple',label[0],label[1])))
  assert keys==channels(split,L)
  bits=cut['coefficientDyadicBits'];assert type(bits)==int and bits>0
  assert len(keys)==len(cut['coefficientNumerators'])==len(cut['columnScaleExponents'])
  gamma=[]
  for a,b,e,num in zip((k[0] for k in keys),(k[1] for k in keys),cut['columnScaleExponents'],cut['coefficientNumerators']):
   assert type(e)==int and e>=0;assert normcouple(a,b,L)>0;gamma.append(Q(int(num),1<<(bits+e)))
  weight=Q(int(cut['weightNumerator']),1<<cut['weightDyadicBits']);assert weight>0
  key=(split,L)
  if key not in QRAW:QRAW[key]=[[Q(0) for _ in keys] for _ in keys];LAY[key]=keys
  for i,a in enumerate(gamma):
   for j in range(i,len(gamma)):
    value=weight*a*gamma[j];QRAW[key][i][j]+=value
    if i!=j:QRAW[key][j][i]+=value
  splitCounts[split]+=1
 assert len(cutIds)==81
 LOOK={key:{ab:i for i,ab in enumerate(keys)} for key,keys in LAY.items()}
 blocks=source['inputBlocks'];variables=source['variables'];target=source['pptFunctional'];targetDen=int(target['denominator']);targetNum=list(map(int,target['numerators']));assert len(targetNum)==len(variables)==4955
 expectedVariables=[[ib,a,b] for ib,B in enumerate(blocks) for a in range(B['multiplicity']) for b in range(a,B['multiplicity'])];assert variables==expectedVariables
 projectionChecks=0;coordinates=0;outputSpinSets={s:sorted(L for ss,L in LOOK if ss==s) for s in ['13','22']}
 position=0
 for ib,B in enumerate(blocks):
  J=B['J2'];p=B['multiplicity'];H=[[int(v) for v in row] for row in B['highestWeight']];ws=BY[4][J];ix={w:i for i,w in enumerate(ws)}
  matrices=[]
  for split in ['13','22']:
   ch=channels(split,J);over=[];nc=[];physicalNorms=[]
   for a,b in ch:
    Z=symcouple(a,b,J);v=[0]*p
    for w,value in Z.items():
     hh=H[ix[w]];factor=value*gnum(w)
     for k,h in enumerate(hh):v[k]+=factor*h
    over.append(v);nc.append(cg_highest(a[1],b[1],J)[3]);physicalNorms.append(normcouple(a,b,J))
   # Exact completeness/normalization of each source embedding in the product basis.
   for a in range(p):
    for b in range(a,p):
     projected=sum(Q(v[a]*v[b],576)/norm for v,norm in zip(over,physicalNorms))
     physical=sum(Q(row[a]*row[b]*gnum(w),24) for w,row in zip(ws,H))
     assert projected==physical;projectionChecks+=1
   m=len(ch);G=[[Q(0) for _ in ch] for _ in ch]
   for u,(a,b) in enumerate(ch):
    for v in range(u,m):
     c,d=ch[v];value=Q(0)
     for key,index in LOOK.items():
      if key[0]!=split:continue
      alpha=(c,b);beta=(a,d)
      if alpha not in index or beta not in index:continue
      w=QRAW[key][index[alpha]][index[beta]]
      if not w:continue
      L=key[1];raw=sixj_raw(a[1],b[1],J,c[1],d[1],L)
      if raw:value+=w*raw*(L+1)*(J+1)
     value/=576*nc[u]*nc[v];G[u][v]=G[v][u]=value
   den=lcm(*(a.denominator for row in G for a in row));GI=[[a.numerator*(den//a.denominator) for a in row] for row in G]
   # Integer matrix products; no FLINT, NumPy or saved numerical overlaps.
   cols=list(zip(*over));middle=[[sum(a*b for a,b in zip(row,col)) for col in cols] for row in GI]
   nums=[[sum(a*middle[k][b] for k,a in enumerate(col))*int(B['denominatorScale']) for b in range(p)] for col in cols]
   matrices.append((nums,den))
  den=lcm(*(d for M,d in matrices));M=[[sum(mat[a][b]*(den//d) for mat,d in matrices) for b in range(p)] for a in range(p)]
  for a in range(p):
   for b in range(a,p):
    assert M[a][b]*(1 if a==b else 2)*targetDen==targetNum[position]*den,(ib,J,a,b)
    position+=1;coordinates+=1
  sixj_raw.cache_clear();print('verifiedPPTsourceblock',ib,J,'coordinates',coordinates,'seconds',time.time()-start,flush=True)
 assert position==4955
 report={'passed':True,'n':25,'copies':4,'splits':[[1,3],[2,2]],'sourceCertificateSHA256':hashlib.sha256(sourceRaw).hexdigest(),'rationalCutsSHA256':hashlib.sha256(cutsRaw).hexdigest(),'positiveDyadicRankOneCuts':len(cutIds),'cutCountsBySplit':dict(splitCounts),'outputDoubledSpins':outputSpinSets,'parentHighestColumnsReconstructed':parentChecks,'parentDescendantNormsChecked':descChecks,'sourceEmbeddingGramEntriesChecked':projectionChecks,'rawSourceFunctionalCoordinatesVerified':coordinates,'method':'Physical integer parent highest/descendant reconstruction; rational CG coupled bases; exact raw6j recoupling; highest-word symmetrization overlaps; exact integer congruences; all saved PPT-functional coefficients independently reconstructed. Standard library only.','scope':'The saved rational functional is a positive sum of 1|3 and2|2 partial-time-reversal compressions for every separable four-copy moment. Combining this with the independently verified source identity proves exclusion.','elapsedSeconds':time.time()-start}
 return report
if __name__=='__main__':
 root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent;r=verify(root);(root/'ppt-functional-independent.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
