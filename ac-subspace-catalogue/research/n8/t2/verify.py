#!/usr/bin/env python3
"""Independent rational certificate for kmax(8,2)=kmax(8,3)=1.
Recoupling is evaluated from rational highest-weight sums, independently of
Wigner 6j / Clebsch-Gordan libraries and numerical optimization.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent
N=8;j=4;d=9

def C(n,k):return comb(n,k) if 0<=k<=n else 0

def h(n,defect,t):
 return sum((F((-1)**q*factorial(n-defect+q)*factorial(n-q),factorial(n-defect)*factorial(n))*C(t,q)*C(n-t,defect-q)for q in range(defect+1)),F(0))/C(n,defect)

def multipole_eigen(L,J):
 return F((2*L+1)*C(N,L),(L+1)*C(N+L+1,N))*sum(((-1)**(t+L)*C(L,t)*C(t+L,t)*h(N,N-J,t)for t in range(L+1)),F(0))

def verify():
 raw=(ROOT/'certificate.json').read_bytes();cert=json.loads(raw);channels=list(range(0,9,2))
 assert cert['copySpinChannels']==channels
 A=[[multipole_eigen(L,J)for J in channels]for L in range(9)]
 assert [[str(x)for x in row]for row in A]==cert['multipoleMatrix']
 beta=[F(2)]+[F(40,3*L*(L+1)-160)for L in range(1,9)]
 assert list(map(str,beta))==cert['beta']
 w={key:F(value)for key,value in cert['identity'].items()};checks=[]
 for col,J in enumerate(channels):
  X=sum(beta[L]*A[L][col]for L in range(9));Y=sum(beta[L]**2*A[L][col]for L in range(9))
  lhs=w['a5']*A[5][col]
  rhs=w['constant']+w['a1']*A[1][col]+w['a2']*A[2][col]+w['tracePQ']*X+w['tracePSquared']*Y
  assert lhs==rhs,(J,lhs,rhs)
  checks.append({'J':J,'exactlyEqual':True,'coefficient':str(lhs)})
 # a1=a2=0, TrPQ=0, TrP²=2 reduces the identity to a positive contradiction.
 assert w['constant']+2*w['tracePSquared']==-96*5831
 assert w['a5']==96*217971
 # Exact Casimir eigenvalues, including every highest-weight tensor rank.
 Jz=s.diag(*[j-r for r in range(d)]);Jp=s.zeros(d)
 for r in range(1,d):Jp[r-1,r]=s.sqrt(r*(N-r+1))
 Jm=Jp.T
 assert Jp*Jm-Jm*Jp==2*Jz
 assert Jz*Jz+(Jp*Jm+Jm*Jp)/2==20*s.eye(d)
 for L in range(9):
  T=Jp**L;factor=20-s.Rational(L*(L+1),2)
  delta=Jz*T*Jz+(Jp*T*Jm+Jm*T*Jp)/2-factor*T
  assert all(s.simplify(z)==0 for z in delta)
  if L:
   qfactor=-1-factor/s.Rational(20,3)
   assert qfactor!=0 and s.Rational(beta[L].numerator,beta[L].denominator)*qfactor==1
 # A pure order-three state supplies a one-dimensional subspace at both orders.
 vector=s.Matrix([s.sympify(z)for z in cert['lowerWitness']['exactCoefficients']]);assert s.simplify((vector.H*vector)[0])==1
 count=0
 for t in[2,3]:
  M=s.Matrix(t+1,N-t+1,lambda a,b:vector[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(N-t,b)/s.binomial(N,a+b)))
  delta=M*M.H-s.eye(t+1)/(t+1)
  assert delta.applyfunc(s.simplify)==s.zeros(t+1);count+=(t+1)**2
 return {'passed':True,'certificateSha256':hashlib.sha256(raw).hexdigest(),'result':{'n':8,'maximumDimensionAtOrder2':1,'maximumDimensionAtOrder3':1,'dimension2Order2Excluded':True},'recouplingMethod':'Rational highest-weight sums, without Wigner/CG routines.','rationalRecouplingEntriesChecked':45,'exactChannelIdentities':checks,'casimirHighestWeightsChecked':9,'lowerWitnessReducedEntriesChecked':count,'contradictionLatex':cert['contradictionLatex'],'scope':cert['scope']}

if __name__=='__main__':
 report=verify();(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
