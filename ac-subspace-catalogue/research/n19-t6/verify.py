from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial,prod
from itertools import combinations_with_replacement
import json,hashlib,flint,time
start=time.time();root=Path(__file__).resolve().parent;raw=(root/'certificate.json').read_bytes();d=json.loads(raw);n=d['n'];r=d['copies'];by={}
for w in combinations_with_replacement(range(n+1),r):by.setdefault(2*sum(w)-n*r,[]).append(w)
A=[[int(x)for x in row]for row in d['selectedRows']];K=[[int(x)for x in row]for row in d['kernel']];prime=1000003
assert flint.nmod_mat(A,prime).rank()==88 and flint.nmod_mat(K,prime).rank()==6
assert flint.fmpz_mat(A)*flint.fmpz_mat(K)==flint.fmpz_mat(88,6)
P=[]
for ib,B in enumerate(d['blocks']):
 J=B['J2'];p=B['multiplicity'];H=[[int(x)for x in row]for row in B['highestWeight']];met=[F(prod(comb(n,i)**w.count(i)*factorial(w.count(i))for i in range(n+1)),factorial(r))for w in by[J]]
 G=[[sum(row[a]*row[b]*g for row,g in zip(H,met))for b in range(p)]for a in range(p)];assert all(G[a][b]==0 for a in range(p)for b in range(p)if a!=b)
 forms=[F(0)]*6
 for v,(bb,a,b)in enumerate(d['variables']):
  if bb!=ib:continue
  c=int(B['denominatorScale'])*(J+1)*G[a][b]*(1 if a==b else 2)
  forms=[x+c*y for x,y in zip(forms,K[v])]
 P.append(forms)
cert=d['probabilityCertificate'];ids=cert['blockIndices'];T=[sum(row[j]for row in P)for j in range(6)];assert T==list(map(F,cert['traceKernelRow']));assert [P[i]for i in ids]==[[F(x)for x in row]for row in cert['kernelRows']];weights=[int(cert['traceWeight'])]+list(map(int,cert['positiveWeights']));assert all(w>0 for w in weights)
assert all(weights[0]*T[j]+sum(w*P[i][j]for i,w in zip(ids,weights[1:]))==0 for j in range(6))
report={'passed':True,'n':19,'t':6,'copies':3,'momentDimension':1540,'blocks':27,'realInvariantVariables':94,'constraintRank':88,'kernelDimension':6,'probabilitySpins2':cert['spins2'],'traceWeight':str(weights[0]),'positiveProbabilityWeights':list(map(str,weights[1:])),'exactKernelIdentityVerified':True,'method':'Discovery-side exact selected row/kernel ranks and rational physical trace/probability identity; independent physical reconstruction is supplied separately.','certificateSHA256':hashlib.sha256(raw).hexdigest(),'elapsedSeconds':time.time()-start}
(root/'exact-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
