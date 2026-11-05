#!/usr/bin/env python3
"""Verify exported exact first-order AC subspaces by compressed spin matrices.
Independent of the constructor's reduced cross-density-matrix implementation.
"""
from pathlib import Path
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parent

def simp(x):return s.simplify(s.expand_complex(x))

def main():
 path=ROOT/'t1-bases.json';data=json.loads(path.read_text());out=[]
 for record in data['subspaces']:
  n,k=record['n'],record['dimension'];j=s.Rational(n,2)
  V=s.Matrix.hstack(*[s.Matrix(list(map(s.sympify,row['exactCoefficients'])))for row in record['basis']])
  assert V.shape==(n+1,k) and (V.H*V).applyfunc(simp)==s.eye(k)
  Jz=s.diag(*[j-r for r in range(n+1)]);Jminus=s.zeros(n+1)
  for r in range(n):Jminus[r+1,r]=s.sqrt((n-r)*(r+1))
  Jplus=Jminus.T
  for operator in [Jz,Jplus,Jminus]:assert (V.H*operator*V).applyfunc(simp)==s.zeros(k),(record['id'],operator)
  assert k==n//3+1
  out.append({'id':record['id'],'n':n,'dimension':k,'orthonormal':True,'compressedSpinOperatorsExactlyZero':True,'entriesChecked':3*k*k})
  print('Verified',record['id'],flush=True)
 # Check the phase-feasibility and scalar-cancellation identities independently.
 X=s.Symbol('X',positive=True);Y=s.Symbol('Y',nonnegative=True);general=[]
 for which in [3,4]:
  if which==3:
   a=8*X/(9*(4*X-7));b=7*(4*X-15)/(18*(4*X-7));c=7*(4*X-3)/(18*(4*X-7));minimum=s.Rational(63,4)
   cosine2=(X+s.Rational(1,4))*b*b/(4*(X-s.Rational(35,4))*a*c)
   assertions=[a+b+c-1,7*a+b-5*c,b*(X-s.Rational(3,4))-c*(X-s.Rational(15,4))]
  else:
   a=16*X/(52*X-135);b=27*(4*X-35)/(8*(52*X-135));c=45*(4*X-3)/(8*(52*X-135));minimum=s.Rational(99,4)
   cosine2=(X+s.Rational(1,4))*b*b/(4*(X-s.Rational(63,4))*a*c)
   assertions=[a+b+c-1,9*a+b-7*c,5*b*(X-s.Rational(3,4))-3*c*(X-s.Rational(35,4))]
  assert all(s.cancel(x)==0 for x in assertions)
  for expression in [a,b,c,1-cosine2]:
   numerator,denominator=s.fraction(s.factor(expression))
   for poly in [numerator,denominator]:
    assert all(coefficient>0 for coefficient in s.Poly(s.expand(poly.subs(X,minimum+Y)),Y).all_coeffs())
  general.append({'centralDimension':which,'normalizationAndMomentIdentities':True,'positiveProbabilitiesAndSineSquaredForEntireDomain':True,'minimumX':str(minimum)})
 result={'allChecksPassed':True,'inputSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'verificationMethod':'Exact compression of Jz,J+,J- in every exported basis, independent of the reduced-matrix constructor.',
         'subspaces':out,'generalCentralBlockChecks':general,'generalExistenceConclusion':'For every n>=4 (and also n=2), there exists a first-order anticoherent subspace of dimension floor(n/3)+1. This alone does not prove global maximality.'}
 (ROOT/'t1-verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print('All',len(out),'exact bases and general phase-feasibility identities verified.')
if __name__=='__main__':main()
