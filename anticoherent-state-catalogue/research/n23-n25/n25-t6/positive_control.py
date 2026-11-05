"""Exact 1|3 and2|2 partial-time-reversal embedding controls for N25, standard library.
No AC property is assumed. This generic complex vector tests the convention.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from collections import Counter
from math import factorial,comb,prod
import json,time
P=Path(__file__).resolve().parent
ZERO=(F(0),F(0))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(x,a):return x[0]*a,x[1]*a
def conj(x):return x[0],-x[1]
def abssq(x):return x[0]*x[0]+x[1]*x[1]
def product(xs):
 out=(F(1),F(0))
 for x in xs:out=mul(out,x)
 return out
def multinomial(w):return factorial(len(w))//prod(factorial(c) for c in Counter(w).values())
def metric(w):return F(prod(comb(25,i) for i in w)*prod(factorial(c) for c in Counter(w).values()),factorial(len(w)))
def reflected(a):return tuple(sorted(25-i for i in a))
def mu(a,b):
 ca,cb=Counter(a),Counter(b);return prod(comb(ca[i]+cb[i],ca[i]) for i in set(a+b))
def verify():
 start=time.time();amp={0:(F(2),F(1)),3:(F(-3),F(2)),11:(F(1),F(-4)),20:(F(2),F(-3)),25:(F(-1),F(-2))};support=set(amp)
 norm=sum(comb(25,i)*abssq(v) for i,v in amp.items())
 words4=list(combinations_with_replacement(sorted(support),4));v4={w:scale(product(amp[i] for i in w),multinomial(w)) for w in words4}
 trace_x=sum(metric(w)*abssq(v) for w,v in v4.items())/norm**4;assert trace_x==1
 controls=[]
 for first in [1,2]:
  second=4-first;choose=comb(4,first)
  va={w:scale(product(amp[i] for i in w),multinomial(w)) for w in combinations_with_replacement(sorted(support),first)}
  vb={w:scale(product(amp[i] for i in w),multinomial(w)) for w in combinations_with_replacement(sorted(support),second)}
  phi={(a,b):scale(mul(conj(va[reflected(a)]),vb[b]),(-1)**sum(a)) for a in sorted(reflected(w) for w in va) for b in vb}
  count=0
  for (a,b),aa in phi.items():
   for (c,d),bb in phi.items():
    cb=reflected(c);ab=reflected(a);w=tuple(sorted(cb+b));v=tuple(sorted(ab+d))
    mapped=scale(mul(v4[w],conj(v4[v])),F((-1)**(sum(a)+sum(c))*mu(cb,b)*mu(ab,d),choose**2*norm**4))
    direct=scale(mul(aa,conj(bb)),F(1,norm**4));assert mapped==direct;count+=1
  rows=0
  for a in combinations_with_replacement(range(26),first):
   for b in combinations_with_replacement(range(26),second):
    assert ((a,b) in phi)==(set(reflected(a))<=support and set(b)<=support);rows+=1
  trace_y=sum(metric(a)*metric(b)*abssq(v) for (a,b),v in phi.items())/norm**4;assert trace_y==1
  controls.append({'split':[first,second],'embeddingDenominator':choose,'densityKernelDenominator':choose**2,'outputProductDimension':rows,'outputNonzeroCoordinates':len(phi),'allNonzeroEntriesCompared':count,'allOutputRowsSupportChecked':rows,'offSupportZeroEquivalenceProved':True,'inputTraceExact':'1','outputTraceExact':'1','outputPositiveByExactRankOneFactorization':True})
 return {'passed':True,'n':25,'copies':4,'splits':[[1,3],[2,2]],'controls':controls,'method':'Generic Gaussian-rational pure state in the binomial basis. Every nonzero partial-time-reversal image entry is compared with direct product-state conjugation; support equivalence proves all other entries vanish.','stateBinomialBasisCoefficients':{str(i):list(map(str,v)) for i,v in amp.items()},'sourceNonzeroOccupationWords':len(v4),'scope':'Convention and positivity control only; this non-anticoherent generic state is not an AC existence witness.','seconds':time.time()-start}
if __name__=='__main__':
 r=verify();(P/'positive-control.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
