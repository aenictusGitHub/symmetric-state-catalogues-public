#!/usr/bin/env python3
"""Exact auxiliary-spin-three AC4 bound; independent physical-matrix audit.
Standard library only. No numerical eigenvalues or discovery code are used.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib
HERE=Path(__file__).resolve().parent

def poly_mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for k,y in enumerate(b):c[i+k]+=x*y
 return c

def poly_add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 while len(c)>1 and c[-1]==0:c.pop()
 return c

def charpoly_tridiagonal(A):
 p=[F(1)];q=[-A[0][0],F(1)]
 for k in range(1,len(A)):
  q,p=poly_add(poly_mul([-A[k][k],F(1)],q),[-A[k][k-1]*A[k-1][k]*c for c in p]),q
 return q

def mm(A,B):
 return [[sum(A[i][k]*B[k][h]for k in range(len(A)))for h in range(len(A))]for i in range(len(A))]

def eye(d):return [[int(i==k)for k in range(d)]for i in range(d)]

def general_identity(n):
 assert n>=16
 j=F(n,2)
 xs=[-3*j-3,-2*j-5,-j-6,-6,j-5,2*j-3,3*j]
 ranks=[j-3+p for p in range(7)];dims=[int(2*L+1)for L in ranks]
 vals=[(x+2*j+5)*(x+6)*(x-j+5)*(x-3*j)for x in xs]
 mu=2*j*(4*j**3-36*j*j+47*j-75)
 h=j-8;assert mu==2*j*(4*h**3+60*h*h+239*h+45)>0
 assert sum(w*f for w,f in zip(dims,vals))==7*(n+1)*mu
 assert vals[2]==-2*j*(j-1)*(2*j+1)*(2*j+3)<0
 assert vals[5]==-2*(j+2)*(j+3)*(2*j+1)*(2*j+3)<0
 assert all(vals[p]==0 for p in [1,3,4,6])
 eig=[v-mu for v in vals]
 assert eig[0]>0 and all(e<0 for e in eig[1:])
 assert dims[0]==n-5 and sum(dims[1:])==6*(n+2)
 return {'n':n,'mean':str(mu),'eigenvalues':[str(x)for x in eig],'multiplicities':dims,'positiveIndex':n-5,'negativeIndex':6*(n+2),'upper':(n-5)//7}

def verify(directory=None):
 root=Path(directory or HERE);raw=(root/'certificate.json').read_bytes();cert=json.loads(raw)
 assert cert['schema']=='spin3-quartic-ac4-inertia-v1' and cert['auxiliarySpin']==3 and cert['minimumN']==16 and cert['minimumACOrder']==4
 assert cert['generalBound']=='7 K(N,t) <= N-5'
 assert cert['generalPolynomialFactors']==['X+2*j+5','X+6','X-j+5','X-3*j']
 assert cert['generalCenteringCoefficientsAscending']==[0,-150,94,-72,8]
 assert cert['positiveCenteringShift8CoefficientsAscending']==[45,239,60,4]
 # Prove the polynomial identity for all j, by exact coefficient arithmetic.
 # The total trace numerator has degree five. Six evaluations certify it
 # identically after clearing its known denominator 7(2j+1).
 for n in [16,17,18,19,20,21]:general_identity(n)
 # Independently build the physical spin9 tensor spin3 matrices in the
 # binomially rescaled magnetic basis, where D=2 J.S is integral.
 n=18;s2=6;blocks={}
 for i in range(n+1):
  for a in range(s2+1):blocks.setdefault(i+a-12,[]).append((i,a))
 tracepowers=[0]*5;partial=[[F(0)for b in range(7)]for a in range(7)]
 charchecks=0;metricchecks=0;counts={L:0 for L in range(6,13)}
 for M,words in sorted(blocks.items()):
  lookup={w:k for k,w in enumerate(words)};d=len(words);D=[[0]*d for _ in range(d)]
  for col,(i,a)in enumerate(words):
   D[col][col]=2*(i-9)*(a-3)
   if i<18 and a>0:D[lookup[i+1,a-1]][col]=(i+1)*(7-a)
   if i>0 and a<6:D[lookup[i-1,a+1]][col]=(19-i)*(a+1)
  metric=[comb(18,i)*comb(6,a)for i,a in words]
  for r in range(d):
   for c in range(d):assert metric[r]*D[r][c]==metric[c]*D[c][r];metricchecks+=1
  expected=[F(1)]
  for L in range(max(6,abs(M)),13):
   expected=poly_mul(expected,[-F(L*(L+1)-102),F(1)]);counts[L]+=1
  assert charpoly_tridiagonal(D)==expected;charchecks+=1
  powers=[eye(d)]
  for k in range(1,5):powers.append(mm(powers[-1],D))
  for k in range(5):tracepowers[k]+=sum(powers[k][a][a]for a in range(d))
  # 16W=D^4-4D^3-2612D^2-9168D+138240 I.
  W=[[F(powers[4][a][b]-4*powers[3][a][b]-2612*powers[2][a][b]-9168*powers[1][a][b]+138240*powers[0][a][b],16)for b in range(d)]for a in range(d)]
  for a,(i,auxa)in enumerate(words):
   for b,(h,auxb)in enumerate(words):
    if i==h:partial[auxa][auxb]+=W[a][b]
 assert counts=={L:2*L+1 for L in range(6,13)}
 assert all(x==0 for row in partial for x in row)
 normalizedmoments=[F(v,133*2**k)for k,v in enumerate(tracepowers)]
 assert normalizedmoments==[1,0,360,-180,226080]
 result=general_identity(18)
 case=cert['n18'];assert case['n']==18 and case['spin']==9 and case['dimension']==133
 assert case['totalSpins']==list(range(6,13)) and case['multiplicities']==result['multiplicities']
 assert case['dotProductEigenvalues']==[-30,-23,-15,-6,4,15,27] and case['quarticRoots']==[-23,-6,4,27]
 assert case['centeringConstant']==6264 and case['polynomialCoefficientsAscending']==[8640,-1146,-653,-2,1]
 assert case['witnessEigenvalues']==list(map(int,result['eigenvalues']))
 assert (case['positiveIndex'],case['negativeIndex'],case['maximumDimensionUpper'])==(13,120,1)
 assert result['eigenvalues']==list(map(str,[319320,-6264,-63720,-6264,-6264,-111600,-6264]))
 report={'passed':True,'certificateSHA256':hashlib.sha256(raw).hexdigest(),'n':18,'minimumACOrder':4,'maximumDimensionUpper':1,'physicalHilbertDimension':133,'physicalMagneticBlocksChecked':charchecks,'weightedHermiticityEntriesChecked':metricchecks,'allCharacteristicPolynomialsExact':True,'partialTraceOverSpinNineExactlyZero':True,'partialTraceEntriesChecked':49,'normalizedPhysicalTraceMoments':[str(x)for x in normalizedmoments],'polynomial':'X^4 - 2 X^3 - 653 X^2 - 1146 X + 8640 I','spectrum':result,'generalBound':{'formula':'7 K(N,t) <= N-5','minimumN':16,'minimumT':4,'domainAllIntegers':True,'positiveMeanExpansion':'2j [4(j-8)^3+60(j-8)^2+239(j-8)+45]','generalTraceIdentityVerifiedByExactPolynomialInterpolation':True},'method':'Independent integral physical spin matrices; weighted Hermiticity, every magnetic-block characteristic polynomial and all partial-trace entries checked exactly. Inertia argument and tensor-rank compression given in proof.tex. Standard library only; no numerical eigensolver.','scope':'Upper bound for arbitrary complex AC subspaces, with no support or symmetry restriction.'}
 return report
if __name__=='__main__':
 report=verify();(HERE/'upper-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
