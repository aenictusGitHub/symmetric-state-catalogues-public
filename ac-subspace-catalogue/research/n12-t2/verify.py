#!/usr/bin/env python3
"""Exact independent proof certificate: N=12, t=2 has maximum dimension2.
Requires SymPy only for the existing lower-bound witness and spin generators.
All upper-bound recoupling and coefficient identities use rational arithmetic.
No numerical optimization or Wigner/Clebsch–Gordan routines are used.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent
N=12;d=13;j=6

def C(n,k):return comb(n,k) if 0<=k<=n else 0

def h(n,defect,t):
 return sum((F((-1)**q*factorial(n-defect+q)*factorial(n-q),factorial(n-defect)*factorial(n))*C(t,q)*C(n-t,defect-q) for q in range(defect+1)),F(0))/C(n,defect)

def multipole_eigen(L,J):
 return F((2*L+1)*C(N,L),(L+1)*C(N+L+1,N))*sum(((-1)**(t+L)*C(L,t)*C(t+L,t)*h(N,N-J,t) for t in range(L+1)),F(0))

raw=(ROOT/'certificate.json').read_bytes();cert=json.loads(raw);channels=list(range(0,13,2))
A=[[multipole_eigen(L,J) for J in channels]for L in range(13)]
assert [[str(x)for x in row]for row in A]==cert['multipoleMatrix']
beta=[F(3)]+[F(28,L*(L+1)-112)for L in range(1,13)]
assert list(map(str,beta))==cert['beta']
c=cert['identity'];a_weights=[F(c[f'a{L}'])for L in [1,2,3]];p_weights=list(map(F,c['pJ']));constant=F(c['constant'])
checks=[]
for column,J in enumerate(channels):
 lhs=-sum(beta[L]*A[L][column]for L in range(13))
 rhs=constant+sum(a_weights[L-1]*A[L][column]for L in [1,2,3])+p_weights[column]
 assert lhs==rhs,(J,lhs,rhs)
 checks.append({'J':J,'left':str(lhs),'right':str(rhs),'exactlyEqual':True})
scale=F(1037300,9)
assert scale*constant==9015
assert scale*a_weights[2]==1031316
assert [scale*w for w in p_weights]==[3936270,20020,0,728365,0,0,0]
assert constant>0 and a_weights[2]>0 and all(w>=0 for w in p_weights)
# Exact spin-6 algebra verifies the Casimir superoperator eigenvalues.
Jz=s.diag(*[j-r for r in range(d)]);Jp=s.zeros(d)
for r in range(1,d):Jp[r-1,r]=s.sqrt(r*(N-r+1))
Jm=Jp.T
assert Jp*Jm-Jm*Jp==2*Jz
assert Jz*Jp-Jp*Jz==Jp
assert Jz*Jz+(Jp*Jm+Jm*Jp)/2==42*s.eye(d)
for L in range(13):
 T=Jp**L;factor=s.Rational(84-L*(L+1),2)
 delta=Jz*T*Jz+(Jp*T*Jm+Jm*T*Jp)/2-factor*T
 assert all(s.simplify(x)==0 for x in delta)
# Existing exact two-dimensional second-order construction, expressed via phases.
v=(7*s.sqrt(285)+2*s.I*s.sqrt(269515))/1045;w=-s.conjugate(v)
u1=s.zeros(d,1);u2=s.zeros(d,1)
for r,z in {3:s.sqrt(38)/8,7:s.Rational(1,4),11:s.sqrt(22)*v/8}.items():u1[r]=z
for r,z in {1:s.sqrt(22)/8,5:w/4,9:s.sqrt(38)*s.conjugate(v)/8}.items():u2[r]=z
B=s.Matrix.hstack(u1,u2)
assert (B.H*B-s.eye(2)).applyfunc(s.simplify)==s.zeros(2)
matrices=[s.Matrix(3,11,lambda a,b:u[a+b]*s.sqrt(s.binomial(2,a)*s.binomial(10,b)/s.binomial(12,a+b)))for u in [u1,u2]]
lower_checks=0
for a,Ma in enumerate(matrices):
 for b,Mb in enumerate(matrices):
  delta=Ma*Mb.H-(s.eye(3)/3 if a==b else s.zeros(3))
  assert delta.applyfunc(lambda z:s.simplify(s.expand_complex(z)))==s.zeros(3)
  lower_checks+=9
report={'passed':True,'result':{'n':12,'order':2,'maximumDimension':2,'dimension3Excluded':True,'dimension2ExactlyAttained':True},'certificateSha256':hashlib.sha256(raw).hexdigest(),'rationalRecouplingEntriesChecked':91,'channelIdentities':checks,'casimirHighestWeightsChecked':13,'lowerBoundExactCrossEntriesChecked':lower_checks,'contradictionLatex':r'0=9015+1031316a_3+3936270p_0+20020p_2+728365p_6>0','nonnegativeQuantities':['a_3: squared Hilbert–Schmidt norm of the rank3 multipole component of a pure complement state','p_J: total-spin-J probabilities of two identical copies of that state'],'scope':'The upper bound applies to every complex pure-state subspace; no support, reality, Majorana multiplicity, symmetry or numerical-search assumption is used.'}
(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
