#!/usr/bin/env python3
"""Exact radical N=17, t=4 witnesses; no numerical search enters the proof."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
SUPPORT=(1,4,7,10,13,16)
def coefficients(sign=1):
 q=s.sqrt(20995)
 v=s.sqrt(12780*q-575935)
 a=(710-9*q-sign*s.I*v)/(756*s.sqrt(15))
 b=2*s.sqrt(21)/27
 c=2*(9*q-2840-sign*2*s.I*v)/(189*s.sqrt(4290))
 state=s.zeros(18,1)
 for r,value in zip(SUPPORT,[a,b,c,c,b,a]):state[r]=value
 return state

def zero(value,label):
 value=s.simplify(s.expand_complex(value))
 assert value==0,(label,str(value))

def split(state,t):
 return s.Matrix(t+1,18-t,lambda a,b:state[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(17-t,b)/s.binomial(17,a+b)))

def verify():
 branches=[]
 for sign in [1,-1]:
  state=coefficients(sign);zero((state.H*state)[0]-1,'normalization')
  weights=[s.Rational(277,1458),s.Rational(28,243),s.Rational(142,729)]
  for r,p in zip(SUPPORT,weights+weights[::-1]):zero(state[r]*s.conjugate(state[r])-p,'probability')
  M=split(state,4);rho=M*M.H
  for a in range(5):
   for b in range(5):zero(rho[a,b]-s.Rational(int(a==b),5),'fourth reduction')
  M=split(state,5);rho5=(M*M.H).applyfunc(lambda z:s.simplify(s.expand_complex(z)))
  expected=s.eye(6)/6
  for a,b,value in [(0,3,-s.sqrt(8398)/9282),(1,4,s.sqrt(20995)/9282),(2,5,-s.sqrt(8398)/9282)]:expected[a,b]=expected[b,a]=value
  for value in rho5-expected:zero(value,'fifth reduction')
  purity=s.simplify(s.trace(rho5*rho5));zero(purity-s.Rational(5443,32487),'fifth purity')
  branches.append({'sign':sign,'normalizedExactly':True,'fourthReduction':'I_5/5','fourthReductionEntriesChecked':25,'fifthReductionEntriesChecked':36,'fifthReductionPurity':str(purity),'purityAboveOneSixth':'19/21658','exactACOrder':4})
 q=s.sqrt(20995);u=710-9*q;v=s.sqrt(12780*q-575935)
 zero(u*u+v*v-1628760,'phase modulus')
 assert (12780*144-575935)>0 and 144**2<20995
 if (HERE/'basis.json').exists():
  source=json.loads((HERE/'basis.json').read_text())
  for expression,value in zip(source['basis'][0]['exactCoefficients'],coefficients()):zero(s.sympify(expression)-value,'exported coefficient')
  basis_hash=hashlib.sha256((HERE/'basis.json').read_bytes()).hexdigest()
 else:basis_hash=None
 report={'passed':True,'n':17,'t':4,'dimension':1,'arithmetic':'Exact rational and radical SymPy expressions; every reduced-matrix entry checked.','branches':branches,'radicandStrictlyPositiveExactly':True,'sourceBasisSha256':basis_hash,'scope':'An analytical construction suggested by the Julia support pattern. Equality or collective-rotation equivalence to the numerical Julia representative is not asserted.'}
 return report
if __name__=='__main__':
 report=verify();(HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
