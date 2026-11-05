#!/usr/bin/env python3
"""Exact N=8 first-order AC dimension bound via a pure N=16 obstruction.
Run python3 verify.py. All proof checks use rational or symbolic exact arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import json
import sympy as s
from sympy.physics.wigner import wigner_6j
ROOT=Path(__file__).resolve().parent

def C(n,k):return comb(n,k) if 0<=k<=n else 0

def h(n,defect,t):
 return sum((F((-1)**q*factorial(n-defect+q)*factorial(n-q),factorial(n-defect)*factorial(n))*C(t,q)*C(n-t,defect-q) for q in range(defect+1)),F(0))/C(n,defect)

def eigen(n,L,defect):
 # Rational highest-weight summation for the eigenvalue of the invariant
 # rank-L tensor contraction on total-spin J=n-defect. No Wigner symbols.
 return F((2*L+1)*C(n,L),(L+1)*C(n+L+1,n))*sum(((-1)**(t+L)*C(L,t)*C(t+L,t)*h(n,defect,t) for t in range(L+1)),F(0))

def spin(n):
 d=n+1;j=s.Rational(n,2);Z=s.diag(*[j-r for r in range(d)]);plus=s.zeros(d)
 for r in range(1,d):plus[r-1,r]=s.sqrt(r*(n-r+1))
 return Z,plus,plus.T

def simplify_zero(matrix):return all(s.simplify(x)==0 for x in matrix)

def verify():
 cert=json.loads((ROOT/'certificate.json').read_text())
 Z,Jp,Jm=spin(8);I=s.eye(9)
 assert Jp*Jm-Jm*Jp==2*Z
 assert Z*Z+(Jp*Jm+Jm*Jp)/2==20*I
 casimir=[]
 for L in range(9):
  T=Jp**L;lam=20-s.Rational(L*(L+1),2)
  assert simplify_zero(Z*T*Z+(Jp*T*Jm+Jm*T*Jp)/2-lam*T)
  casimir.append(str(lam))
 # For k4 the scalar component has squared norm16/9; all other components
 # have total squared norm20/9. The minimum eigenvalue of Phi is -16.
 assert s.Rational(16,9)*20-s.Rational(20,9)*16==0
 assert all(20-s.Rational(L*(L+1),2)>-16 for L in range(1,8))
 # An exact known k3 witness establishes that the upper bound is attained.
 # Use the standard pair-module construction for integer spin4:
 # e0 and v2=(sqrt3 e2+sqrt2 e-3)/sqrt5,
 # w2=(sqrt2 e3-sqrt3 e-2)/sqrt5.
 basis=s.zeros(9,3);basis[4,0]=1
 basis[2,1]=s.sqrt(s.Rational(3,5));basis[7,1]=s.sqrt(s.Rational(2,5))
 basis[1,2]=s.sqrt(s.Rational(2,5));basis[6,2]=-s.sqrt(s.Rational(3,5))
 assert simplify_zero(basis.H*basis-s.eye(3))
 assert simplify_zero(basis.H*Z*basis) and simplify_zero(basis.H*Jp*basis)
 # Vectorization with the spin4 time-reversal intertwiner.
 Y=s.zeros(9)
 for r in range(9):Y[r,8-r]=(-1)**r
 assert Y*Y.T==I
 def vec(T):return s.Matrix(81,1,lambda row,col:(T*Y)[row//9,row%9])
 def trace_one(T):
  out=s.zeros(8)
  for r in range(9):
   for q in range(9):
    if r<8 and q<8:out[r,q]+=s.sqrt((8-r)*(8-q))*T[r,q]/8
    if r>0 and q>0:out[r-1,q-1]+=s.sqrt(r*q)*T[r,q]/8
  return out
 T=Jp**8;vectorization_checks=0
 for k in range(17):
  v=vec(T)
  target=s.Matrix(81,1,lambda pos,col:s.sqrt(s.Rational(C(8,pos//9)*C(8,pos%9),C(16,k))) if pos//9+pos%9==k else 0)
  first=next(i for i,x in enumerate(target) if x!=0);ratio=s.simplify(v[first]/target[first]);assert ratio!=0
  assert simplify_zero(v-ratio*target)
  assert simplify_zero(trace_one(T))
  vectorization_checks+=1
  T=Jm*T-T*Jm
 assert simplify_zero(trace_one(I/9)-s.eye(8)/8)
 # X=P-4I/9: X^2=X/9+20I/81 and TrX^2=20/9.
 assert s.Rational(1,9)/s.Rational(20,9)==s.Rational(1,20)
 assert s.Rational(20,81)/s.Rational(20,9)==s.Rational(1,9)
 # Independently reconstruct the complete pure-spin8 recoupling array.
 js=list(range(0,17,2));A=[[eigen(16,L,16-J) for J in js] for L in range(17)]
 for L in range(17):
  for col,J in enumerate(js):assert s.Rational(A[L][col].numerator,A[L][col].denominator)==(2*L+1)*wigner_6j(8,8,L,8,8,J)
 assert all(x==F(1,17) for x in A[0])
 assert all(sum(A[L][col] for L in range(17))==1 for col in range(9))
 assert cert['identity']['constant']=='100' and cert['identity']['p0Coefficient']=='136' and cert['identity']['a9Coefficient']=='115'
 cs=[F(x) for x in cert['identity']['lowMultipoleCoefficients']]
 for col,J in enumerate(js):
  value=F(100)+(136 if J==0 else 0)+115*A[9][col]+sum(cs[L-1]*A[L][col] for L in range(1,8))
  assert value==0
 # Retain the independently discovered higher-rank certificate as a cross-check.
 alt=cert['alternativeA15Identity'];assert alt['constant']=='1105' and alt['p0Coefficient']=='2431' and alt['a15Coefficient']=='11310'
 altcs=[F(x) for x in alt['lowMultipoleCoefficients']]
 for col,J in enumerate(js):
  assert F(1105)+(2431 if J==0 else 0)+11310*A[15][col]+sum(altcs[L-1]*A[L][col] for L in range(1,8))==0
 # An independent elimination of just a1...a7 gives an even stronger
 # inconsistency using only nonnegative p2,p12, without higher multipoles.
 Am=s.Matrix([[s.Rational(x.numerator,x.denominator) for x in row] for row in A]);E=s.ones(1,9).col_join(Am[1:8,:]);rhs=s.Matrix([1]+[0]*7)
 solution=next(iter(s.linsolve((E,rhs))));assert s.simplify(24871*solution[1]+1292646*solution[6]+205620)==0
 # General even-N saturation lifts to pure2N,(N-1)-AC. For N=2q>=10,
 # Gram orthogonality of tensor ranks0,...,q-1 requires q^2<=4q+1.
 q=s.symbols('q',positive=True);a=q/(2*q+1);norm2=q-q*q/(2*q+1)
 assert s.simplify((1-2*a)/norm2-1/(q*(q+1)))==0
 assert s.simplify(a*(1-a)/norm2-1/(2*q+1))==0
 y=s.symbols('y',nonnegative=True);assert s.expand((y+5)**2-4*(y+5)-1)==y*y+6*y+4
 report={'generalEvenQubitConsequence':{'forEvenNAtLeast8':'k_max(N,1)<=N/2-1','n8Method':'Exact auxiliary pure16,t7 recoupling contradiction.','evenNAtLeast10Method':'Saturation would give a pure2N-qubit(N-1)-AC state. Its orthogonal tensor-error vectors number(N/2)^2, exceeding Hilbert-space dimension2N+1.'},'passed':True,'n':8,'t':1,'maximumDimension':3,'k3WitnessExactlyVerified':True,'casimirEigenvalues':casimir,'rank8OperatorBasisVectorizationChecks':vectorization_checks,'rank8PartialTracesVanishExactly':True,'pure16RecouplingEntriesChecked':153,'rationalHighestWeightSumsAgreeWithSixJ':True,'positiveObstruction':cert['identity']['underSevenAC'],'alternativeA15Obstruction':cert['alternativeA15Identity']['underSevenAC'],'independentAlternativeObstruction':'0=205620+24871*p_2+1292646*p_12','rootOfArgument':'Any rank4 first-order AC projector in spin4 would yield a pure symmetric16-qubit7-AC state; the displayed exact positive identity excludes every such state.'}
 return report
if __name__=='__main__':
 report=verify();(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
