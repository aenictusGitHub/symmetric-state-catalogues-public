#!/usr/bin/env python3
"""Independent exact audit of the wedge/four-copy contraction constraints."""
import sympy as s,json
from pathlib import Path
from itertools import combinations
ROOT=Path(__file__).resolve().parent
n=7;d=8;j=s.Rational(7,2)
u=s.zeros(d,1);v=s.zeros(d,1)
u[0]=s.sqrt(s.Rational(3,10));u[5]=s.sqrt(s.Rational(7,10))
v[2]=s.sqrt(s.Rational(7,10));v[7]=-s.sqrt(s.Rational(3,10))
B=u.row_join(v);assert B.H*B==s.eye(2)
Z=s.diag(*[j-r for r in range(d)]);Jp=s.zeros(d)
for r in range(1,d):Jp[r-1,r]=s.sqrt(r*(n-r+1))
Jm=Jp.T;c=j*(j+1)/3
Ts=[('T10',Z),('T11',Jp),('T1m1',Jm),('T20',Z*Z-c*s.eye(d)),('T21',Z*Jp+Jp*Z),('T2m1',Z*Jm+Jm*Z),('T22',Jp*Jp),('T2m2',Jm*Jm)]
W=(u*v.T-v*u.T)/s.sqrt(2);assert s.simplify(s.trace(W*W.H))==1
w={(a,b):W[a,b] for a in range(d) for b in range(d) if W[a,b]!=0}
psi={ab+cd:s.simplify(a*b) for ab,a in w.items() for cd,b in w.items()}
assert s.simplify(sum(z*s.conjugate(z) for z in psi.values()))==1
for (a,b,c_,d_),z in psi.items():
 assert psi.get((b,a,c_,d_),0)==-z
 assert psi.get((a,b,d_,c_),0)==-z
 assert psi.get((c_,d_,a,b),0)==z
for a,b,c_,d_ in combinations(range(d),4):
 assert s.simplify(W[a,b]*W[c_,d_]-W[a,c_]*W[b,d_]+W[a,d_]*W[b,c_])==0
results=[]
for name,T in Ts:
 assert all(s.simplify(z)==0 for z in B.H*T*B)
 # Direct physical-index definition. No precomputed wedge contraction map.
 contraction=s.Matrix(d,d,lambda b,e:sum(T[a,c_]*W[c_,b]*s.conjugate(W[a,e]) for a in range(d) for c_ in range(d)))
 assert all(s.simplify(z)==0 for z in contraction)
 # Independently sum the nonzero entries of the four-copy pure density.
 output={}
 for (a,b,c_,d_),left in psi.items():
  for (e,f,g,h),right in psi.items():
   value=T[e,a]*left*s.conjugate(right)
   if value!=0:
    key=(b,c_,d_,f,g,h);output[key]=output.get(key,0)+value
 assert all(s.simplify(z)==0 for z in output.values())
 results.append({'tensor':name,'codeCompressionZero':True,'oneWedgeContractionZero':True,'fourCopyFirstLegContractionZero':True,'fourCopyOutputEntriesChecked':len(output)})
# Check the normalization convention used by the canonical weighted-word basis.
weighted_norm=0;canonical_count=0
for (a,b,c_,d_),z in psi.items():
 if not (a<b and c_<d_ and (a,b)<=(c_,d_)):continue
 weight=s.binomial(n,a)*s.binomial(n,b)*s.binomial(n,c_)*s.binomial(n,d_)
 raw=s.simplify(z/s.sqrt(weight));orbit=4 if (a,b)==(c_,d_) else 8
 weighted_norm+=orbit*weight*raw*s.conjugate(raw);canonical_count+=1
assert s.simplify(weighted_norm)==1
report={'passed':True,'witness':'N7,t2,k2: u=sqrt(3/10)D0+sqrt(7/10)D5; v=sqrt(7/10)D2-sqrt(3/10)D7.','wedgeNormalization':'w=(u tensor v-v tensor u)/sqrt2','allRank1And2TensorComponentsChecked':8,'contractionDefinition':'[Tr_1((T tensor I)|w><w|)]_(b,e)=sum_(a,c) T_(a,c) w_(c,b) conjugate(w_(a,e)).','young22SupportVerified':'Antisymmetric in each pair, symmetric under pair exchange, and all Plucker relations vanish.','fourCopyNonzeroCoefficients':len(psi),'canonicalWeightedWordCount':canonical_count,'weightedOrbitNormalizationVerified':True,'components':results,'conceptualIdentity':'If W=(uv^T-vu^T)/sqrt2, the contraction is (uu*/2)<v|T|v>+(vv*/2)<u|T|u>-(uv*/2)<u|T|v>-(vu*/2)<v|T|u>. It vanishes iff the complete2x2 compression B*TB vanishes. The four-copy contraction factors as this matrix tensor |w><w|.'}
(ROOT/'wedge-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
