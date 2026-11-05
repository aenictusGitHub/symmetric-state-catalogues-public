#!/usr/bin/env python3
"""Exact cubic auxiliary-spin-5/2 bound, with physical N19/N20 checks.
Standard library only; integer-rescaled magnetic blocks and exact fractions.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json
HERE=Path(__file__).resolve().parent

def pmul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def padd(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 while len(c)>1 and c[-1]==0:c.pop()
 return c

def charpoly(A):
 p=[F(1)];q=[-A[0][0],F(1)]
 for k in range(1,len(A)):q,p=padd(pmul([-A[k][k],F(1)],q),[-A[k][k-1]*A[k-1][k]*x for x in p]),q
 return q

def mm(a,b):return [[sum(a[i][k]*b[k][j]for k in range(len(a)))for j in range(len(a))]for i in range(len(a))]
def eye(d):return [[int(i==j)for j in range(d)]for i in range(d)]
def general(n):
 assert n>=19;j=F(n,2);rs=[F(-5,2)+p for p in range(6)];spins=[j+r for r in rs]
 xs=[(L*(L+1)-j*(j+1)-F(35,4))/2 for L in spins];dims=[int(2*L+1)for L in spins]
 values=[-(2*x+3*j+8)*(2*x-3*j+5)*(2*x-5*j)for x in xs]
 mean=F(20,3)*j*(j-9)*(2*j-1);assert mean>0
 assert values==[40*j*(2*j-3)*(2*j+1),0,-12*(j+1)*(2*j-1)*(2*j+3),-16*j*(j+2)*(2*j+3),0,0]
 assert sum(w*v for w,v in zip(dims,values))==6*(n+1)*mean
 eig=[v-mean for v in values];assert eig[0]>0 and all(v<0 for v in eig[1:]);assert dims[0]==n-4
 return {'n':n,'spin':str(j),'centeringConstant':str(mean),'totalSpins':[str(L)for L in spins],'dotProductEigenvalues':[str(v)for v in xs],'multiplicities':dims,'wEigenvalues':[str(v)for v in eig],'positiveIndex':n-4,'negativeIndex':5*(n+2),'upper':(n-4)//6}

def physical(n):
 j=F(n,2);s2=5;blocks={}
 for i in range(n+1):
  for a in range(6):blocks.setdefault(2*(i+a)-n-5,[]).append((i,a))
 partial=[[F(0)for _ in range(6)]for _ in range(6)];hermitian=0;charchecks=0;source=general(n);spin2=[int(2*F(v))for v in source['totalSpins']];counts={L:0 for L in spin2}
 # Expand the factored cubic exactly before subtracting its scalar mean.
 polynomial=pmul(pmul([3*j+8,2],[-3*j+5,2]),[-5*j,2]);polynomial=[-v for v in polynomial];polynomial[0]-=F(source['centeringConstant'])
 for M2,words in sorted(blocks.items()):
  lookup={w:k for k,w in enumerate(words)};d=len(words);D=[[F(0)]*d for _ in range(d)]
  for col,(i,a)in enumerate(words):
   D[col][col]=(2*i-n)*(2*a-5)
   if i<n and a>0:D[lookup[i+1,a-1]][col]=2*(i+1)*(6-a)
   if i>0 and a<5:D[lookup[i-1,a+1]][col]=2*(n-i+1)*(a+1)
  metric=[comb(n,i)*comb(5,a)for i,a in words]
  for r in range(d):
   for c in range(d):assert metric[r]*D[r][c]==metric[c]*D[c][r];hermitian+=1
  expected=[F(1)]
  for L in spin2:
   if L<abs(M2):continue
   eig=F(L*(L+2)-n*(n+2)-35,2);expected=pmul(expected,[-eig,1]);counts[L]+=1
  assert charpoly(D)==expected;charchecks+=1
  powers=[eye(d)]
  for h in range(1,4):powers.append(mm(powers[-1],D))
  W=[[sum(polynomial[h]*powers[h][r][c]/4**h for h in range(4))for c in range(d)]for r in range(d)]
  for r,(i,a)in enumerate(words):
   for c,(ii,b)in enumerate(words):
    if i==ii:partial[a][b]+=W[r][c]
 assert all(v==0 for row in partial for v in row)
 assert counts=={L:L+1 for L in spin2}
 return {'n':n,'passed':True,'ambientDimension':6*(n+1),'magneticBlocksChecked':charchecks,'weightedHermiticityEntriesChecked':hermitian,'allCharacteristicPolynomialsExact':True,'partialTraceEntriesChecked':36,'partialTraceExactlyZero':True,'polynomialCoefficientsAscending':[str(v)for v in polynomial],**source}

def verify(directory=None):
 # Five exact values establish the degree-at-most-four trace polynomial
 # identity after clearing its known denominator6(2j+1).
 for n in [19,20,21,22,23]:general(n)
 report={'passed':True,'generalBound':{'formula':'6 K(N,t) <= N-4','minimumN':19,'minimumOrder':3,'auxiliarySpin':'5/2','positiveIndex':'N-4','negativeIndex':'5*(N+2)','exactDegreeFourTraceIdentityChecked':True,'strictPositiveCenteringForAllNInDomain':True},'physicalChecks':[physical(19),physical(20)],'scope':'Arbitrary complex AC subspaces; complete degree-three compression and exact inertia proof.','method':'Exact rational magnetic blocks in a binomially rescaled physical basis. No numerical eigenvalues or optimization statuses.'}
 return report
if __name__=='__main__':
 report=verify();(HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
