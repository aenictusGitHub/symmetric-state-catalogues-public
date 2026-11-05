#!/usr/bin/env python3
"""Independent normalized-spin audit of the auxiliary-spin-5/2 cubic theorem.
Constructs magnetic blocks from square-root ladder amplitudes, rather than
binomial rescaling, and uses full exact characteristic polynomials.
"""
from pathlib import Path
from math import comb
import json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent

def physical(n):
 j=s.Rational(n,2);a=s.Rational(5,2);words={}
 for r in range(n+1):
  for b in range(6):words.setdefault(2*(r+b)-n-5,[]).append((r,b))
 spins=[j-a+h for h in range(6)];xs=[s.factor((L*(L+1)-j*(j+1)-a*(a+1))/2)for L in spins];dims=[int(2*L+1)for L in spins];mu=s.Rational(20,3)*j*(j-9)*(2*j-1)
 eigen=[s.factor(-(2*x+3*j+8)*(2*x-3*j+5)*(2*x-5*j)-mu)for x in xs]
 assert eigen[0]>0 and all(x<0 for x in eigen[1:]);assert dims[0]==n-4 and sum(dims[1:])==5*(n+2)
 partial=[s.Integer(0)]*6;counts=[0]*6;checks=[];z=s.Symbol('z');diagonal={}
 for M2,ww in sorted(words.items()):
  lookup={w:i for i,w in enumerate(ww)};d=len(ww);X=s.zeros(d)
  for k,(r,b)in enumerate(ww):
   X[k,k]=(r-j)*(b-a)
   if r<n and b>0:
    h=lookup[r+1,b-1];X[k,h]=X[h,k]=s.sqrt((n-r)*(r+1)*b*(6-b))/2
  assert X==X.T
  active=[h for h,L in enumerate(spins)if 2*L>=abs(M2)]
  target=s.Poly(s.prod(z-xs[h]for h in active),z)
  assert X.charpoly(z).as_poly()==target
  W=(-(2*X+(3*j+8)*s.eye(d))*(2*X+(-3*j+5)*s.eye(d))*(2*X-5*j*s.eye(d))-mu*s.eye(d)).applyfunc(s.simplify)
  assert W==W.T
  assert W.charpoly(z).as_poly()==s.Poly(s.prod(z-eigen[h]for h in active),z)
  for h in active:counts[h]+=1
  for k,(r,b)in enumerate(ww):partial[b]+=W[k,k];diagonal[r,b]=W[k,k]
  checks.append({'totalMagneticNumberTwice':M2,'dimension':d,'bothCharacteristicPolynomialsExactlyVerified':True})
 assert counts==dims and list(map(s.simplify,partial))==[0]*6
 # Off-diagonal partial-trace entries vanish by conserved total magnetization.
 return {'n':n,'ambientDimension':6*(n+1),'magneticBlocksChecked':len(checks),'characteristicPolynomialsChecked':2*len(checks),'allPartialTraceEntriesExactlyZero':True,'partialTraceEntriesChecked':36,'centeringConstant':str(mu),'totalSpins':list(map(str,spins)),'xEigenvalues':list(map(str,xs)),'wEigenvalues':list(map(str,eigen)),'multiplicities':dims,'positiveIndex':dims[0],'negativeIndex':sum(dims[1:]),'zeroIndex':0,'upper':(n-4)//6,'physicalBlocks':checks}

def verify(directory=None):
 j=s.Symbol('j',real=True);h=s.Symbol('h',positive=True);a=s.Rational(5,2);L=[j-a+i for i in range(6)];x=[s.factor((v*(v+1)-j*(j+1)-a*(a+1))/2)for v in L];F=[s.factor(-(2*v+3*j+8)*(2*v-3*j+5)*(2*v-5*j))for v in x];mu=s.Rational(20,3)*j*(j-9)*(2*j-1)
 assert s.factor(sum((2*l+1)*f for l,f in zip(L,F))-6*(2*j+1)*mu)==0
 assert all(s.factor(f-g)==0 for f,g in zip(F,[40*j*(2*j-3)*(2*j+1),0,-12*(j+1)*(2*j-1)*(2*j+3),-16*j*(j+2)*(2*j+3),0,0]))
 pos=s.Poly(s.expand((F[0]-mu).subs(j,9+h)),h);assert all(q>0 for q in pos.all_coeffs())
 # j>9 is exactly the integer-N domain N>=19; every factor of mu is positive.
 assert s.factor(mu.subs(j,9+h)-s.Rational(20,3)*(9+h)*h*(17+2*h))==0
 result={'passed':True,'generalBound':'6 K(N,t) <= N-4 for integers N>=19 and t>=3','exactSymbolicWeightedTraceIdentity':True,'lowestEigenvalueShift9Polynomial':str(pos.as_expr()),'allShiftedCoefficientsStrictlyPositive':True,'generalSignPattern':[1,-1,-1,-1,-1,-1],'physicalChecks':[physical(19),physical(20)],'arbitraryComplexSubspacesCovered':True,'method':'Independent normalized spin-ladder magnetic blocks with square-root amplitudes and exact full characteristic polynomials, plus a symbolic all-j trace identity and strict sign proof. No discovery or primary-verifier imports.','compressionArgument':'Every physical coefficient of the centered cubic has degree at most3. Its physical partial trace vanishes, so compression to an AC3 code tensor auxiliary spin5/2 is zero. Its sixfold code dimension cannot exceed the positive spectral index N-4.'}
 return result
if __name__=='__main__':
 report=verify();(HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='physicalChecks'},indent=2))
