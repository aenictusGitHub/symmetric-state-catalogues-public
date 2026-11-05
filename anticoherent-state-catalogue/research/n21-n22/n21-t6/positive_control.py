#!/usr/bin/env python3
"""Independent direct symmetric embedding/PT control of the 1|2 map.
No discovery code imported. Product moment identity is proved algebraically;
numerics evaluate every same-weight matrix entry and random cross-weight ones.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
from math import factorial,sqrt
import json,numpy as np,sympy as s
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent

def symmetric_power(c,r):
 out={}
 for w in combinations_with_replacement(range(len(c)),r):
  mult=factorial(r)
  for count in Counter(w).values():mult/=factorial(count)
  out[w]=sqrt(mult)*np.prod([c[i]for i in w])
 return out

def check(c,label):
 n=len(c)-1;r=3;d=n+1
 v4=symmetric_power(c,r);v3=symmetric_power(c,r-1);triples=list(v3);ti={a:i for i,a in enumerate(triples)}
 # E is obtained from symmetrized THREE-copy amplitudes and the physical
 # embedding coefficient, rather than from a product formula.
 E=np.empty((len(triples),d),complex)
 for ia,a in enumerate(triples):
  for i in range(d):E[ia,i]=sqrt((a.count(i)+1)/r)*v4[tuple(sorted(a+(i,)))]
 embed_expected=np.array(list(v3.values()))[:,None]*c[None,:]
 embedding_error=float(np.max(np.abs(E-embed_expected)))
 blocks={}
 for i in range(d):
  for a in triples:blocks.setdefault(2*(i+sum(a))-n*r,[]).append((i,ti[a]))
 theta=(-1.)**np.arange(d)*c[::-1].conj();amp3=np.array(list(v3.values()))
 maxerr=0.;count=0;trace=0.;fullrows=[]
 for M,words in blocks.items():
  arr=np.array(words);i=arr[:,0];a=arr[:,1]
  # Direct first-leg transpose then the antiunitary intertwiner.
  Y=(-1.)**(i[:,None]+i[None,:])*E[a[:,None],n-i[None,:]]*E[a[None,:],n-i[:,None]].conj()
  eta=theta[i]*amp3[a];expected=np.outer(eta,eta.conj())
  maxerr=max(maxerr,float(np.max(np.abs(Y-expected))));count+=Y.size;trace+=np.trace(Y).real
  fullrows.extend(words)
 rng=np.random.default_rng(8391);cross=0.
 for _ in range(10000):
  i,a=fullrows[rng.integers(len(fullrows))];j,b=fullrows[rng.integers(len(fullrows))]
  mapped=(-1.)**(i+j)*E[a,n-j]*E[b,n-i].conjugate()
  expected=theta[i]*amp3[a]*theta[j].conjugate()*amp3[b].conjugate()
  cross=max(cross,float(abs(mapped-expected)))
 assert embedding_error<1e-14 and maxerr<1e-14 and cross<1e-14 and abs(trace-1)<1e-12
 return dict(label=label,n=n,copies=r,symmetricInputDimension=len(v4),outputDimension=d*len(v3),weightBlocks=len(blocks),sameWeightMatrixEntriesChecked=count,crossWeightMatrixEntriesChecked=10000,maximumEmbeddingError=embedding_error,maximumSameWeightMapError=maxerr,maximumCrossWeightMapError=cross,outputTrace=float(trace),outputEqualsPositiveRankOneProjector=True)

rng=np.random.default_rng(48391);reports=[]
for n in [2,21]:
 c=rng.normal(size=n+1)+1j*rng.normal(size=n+1);c/=np.linalg.norm(c)
 reports.append(check(c,'generic complex non-AC state'))

n=21;order=4;exact=[s.S.Zero]*22
for r,a in {0:s.sqrt(915),21:s.sqrt(915),4:13*s.sqrt(14),17:-13*s.sqrt(14),8:17*s.sqrt(7),13:17*s.sqrt(7)}.items():exact[r]=a/(4*s.sqrt(663))
assert s.simplify(sum(z*s.conjugate(z)for z in exact))==1
checks=[]
for q in range(1,5):
 M=s.Matrix(q+1,n-q+1,lambda a,b:exact[a+b]*s.sqrt(s.Rational(s.binomial(q,a)*s.binomial(n-q,b),s.binomial(n,a+b))))
 delta=(M*M.conjugate().T-s.eye(q+1)/(q+1)).applyfunc(s.simplify);assert all(z==0 for z in delta)
 checks.append({'t':q,'maximallyMixedVerifiedExactly':True})
c=np.array([complex(z.evalf(25))for z in exact]);z=np.diag(np.arange(22)-21/2);up=np.zeros((22,22))
for i in range(21):up[i+1,i]=np.sqrt((21-i)*(i+1))
x=(up+up.T)/2;y=(up-up.T)/2j;c=expm(-1j*(.43*x-.79*y+.31*z))@c;reports.append(check(c,'generic rotation of the exact N21 order-four witness'))
report={'passed':True,'exactN21OrderFourChecks':checks,'controls':reports,'algebraicIdentity':'For every pure product moment rho^tensor3, first-leg transpose followed by F|i>=(-1)^i|N-i> gives the positive rank-one projector on Theta(psi) tensor psi^tensor2. The symmetric embedding and all same-weight entries plus random cross-weight entries are independently checked.','scope':'Positive-map convention controls and exact lower-order witness check; this control is not a premise of the exact AC6 exclusion.'}
(HERE/'positive-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
