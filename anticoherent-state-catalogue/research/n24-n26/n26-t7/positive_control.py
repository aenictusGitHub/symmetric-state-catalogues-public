#!/usr/bin/env python3
"""Exact independent physical control of the 1|2 time-reversal embedding.
Standard library only. The test uses a sparse generic complex pure spin-13
vector, without rotational symmetry or anticoherence.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from collections import Counter
from math import comb, factorial, prod
from pathlib import Path
import json

ZERO=(F(0),F(0))
def add(x,y): return x[0]+y[0],x[1]+y[1]
def mul(x,y): return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(x,a): return x[0]*a,x[1]*a
def conj(x): return x[0],-x[1]
def abssq(x): return x[0]*x[0]+x[1]*x[1]
def product(v):
 out=(F(1),F(0))
 for x in v: out=mul(out,x)
 return out

def averaged_word_factor(w):
 return factorial(len(w))//prod(factorial(c) for c in Counter(w).values())
def metric(w):
 return F(prod(comb(26,i) for i in w)*prod(factorial(c) for c in Counter(w).values()),factorial(len(w)))

def verify():
 n=26
 amplitudes={0:(F(2),F(1)),3:(F(-3),F(2)),11:(F(1),F(-4)),20:(F(2),F(-3)),26:(F(-1),F(-2))}
 support=set(amplitudes)
 norm=sum(comb(n,i)*abssq(v) for i,v in amplitudes.items())
 triples=list(combinations_with_replacement(sorted(support),3))
 pairs=list(combinations_with_replacement(sorted(support),2))
 triple_vector={w:scale(product(amplitudes[i] for i in w),averaged_word_factor(w)) for w in triples}
 # Direct tensor product followed by conjugation of the first factor and F.
 phi={(i,a):scale(mul(conj(amplitudes[n-i]),product(amplitudes[j] for j in a)),(-1)**i*averaged_word_factor(a)) for i in sorted(n-j for j in support) for a in pairs}
 checked=0
 for (i,a),va in phi.items():
  for (j,b),vb in phi.items():
   o=tuple(sorted(a+(n-j,)));p=tuple(sorted(b+(n-i,)))
   x=scale(mul(triple_vector[o],conj(triple_vector[p])),F(1,norm**3))
   # Exactly the certificate's 1|2 map, evaluated on a physical complex state.
   mapped=scale(x,F((-1)**(i+j)*(a.count(n-j)+1)*(b.count(n-i)+1),9))
   direct=scale(mul(va,conj(vb)),F(1,norm**3))
   assert mapped==direct,(i,a,j,b,mapped,direct)
   checked+=1
 # Exhaust all product rows to prove the omitted entries vanish on both sides:
 # o,p can lie in the triple support iff the two product rows lie in phi.
 allrows=[(i,a) for i in range(n+1) for a in combinations_with_replacement(range(n+1),2)]
 for row in allrows:
  i,a=row
  assert (row in phi)==(n-i in support and set(a)<=support)
 trace_x=sum(metric(w)*abssq(v) for w,v in triple_vector.items())/norm**3
 trace_y=sum(comb(n,i)*metric(a)*abssq(v) for (i,a),v in phi.items())/norm**3
 assert trace_x==trace_y==1
 assert any(v[1] for v in amplitudes.values())
 return {'passed':True,'n':n,'copies':3,'split':[1,2],
  'control':'Sparse generic complex pure state; exact Gaussian-rational arithmetic; direct physical product-state partial time reversal compared with the certificate embedding.',
  'stateBinomialBasisCoefficients':{str(i):[str(v[0]),str(v[1])] for i,v in amplitudes.items()},
  'squaredNormInPhysicalMetric':str(norm),'inputNonzeroOccupationWords':len(triple_vector),'outputProductDimension':len(allrows),'outputNonzeroProductCoordinates':len(phi),
  'allNonzeroMappedEntriesCompared':checked,'allOutputRowsSupportChecked':len(allrows),'offSupportZeroEquivalenceProved':True,
  'inputTraceExact':str(trace_x),'outputTraceExact':str(trace_y),'outputPositiveByExactRankOneFactorization':True,
  'scope':'Convention and positivity control only; this non-anticoherent state is not an existence witness for any AC order.'}

if __name__=='__main__':
 result=verify();path=Path(__file__).with_name('positive-control.json');path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
