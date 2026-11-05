#!/usr/bin/env python3
"""Independent direct symmetric embedding/PT control of the 1|3 map.
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
 n=len(c)-1;r=4;d=n+1
 v4=symmetric_power(c,4);v3=symmetric_power(c,3);triples=list(v3);ti={a:i for i,a in enumerate(triples)}
 # E is obtained from symmetrized FOUR-copy amplitudes and the physical
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
for n in [2,15]:
 c=rng.normal(size=n+1)+1j*rng.normal(size=n+1);c/=np.linalg.norm(c)
 reports.append(check(c,'generic complex non-AC state'))
cat=HERE.parents[1]/'data/catalogue.json'
D=json.loads(cat.read_text());rec=next(x for x in D['subspaces']if x['id']=='radical-n17-k1-t4')
exact=[s.sympify(x)for x in rec['basis'][0]['exactCoefficients']];n=17;d=18;j=s.Rational(n,2);exactchecks=[]
for q in range(5):
 for p in range(5-q):
  value=sum(s.conjugate(exact[i+q])*exact[i]*(i-j)**p*s.sqrt(s.Rational(factorial(n-i)*factorial(i+q),factorial(n-i-q)*factorial(i)))for i in range(d-q))
  target=sum((i-j)**p for i in range(d))/d if q==0 else s.Integer(0)
  error=s.simplify(s.expand_complex(value-target));assert error==0,(q,p,error)
  exactchecks.append(dict(raisingPower=q,diagonalPower=p,identityVerifiedExactly=True))
c=np.array([complex(x.evalf(25))for x in exact])
z=np.diag(np.arange(d)-n/2);up=np.zeros((d,d))
for i in range(d-1):up[i+1,i]=np.sqrt((n-i)*(i+1))
x=(up+up.T)/2;y=(up-up.T)/2j
c=expm(-1j*(.43*x-.79*y+.31*z))@c
reports.append(check(c,'generic rotation of known exact N17 AC4 state'))
report=dict(passed=True,exactN17AC4MomentChecks=exactchecks,sourceSubspaceId=rec['id'],controls=reports,algebraicIdentity='The normalized symmetric embedding has E[(i,a),o]=sqrt((count_a(i)+1)/4) when o=sort(i,a). Applying first-leg transpose and U|i>=(-1)^i|N-i> gives the stated map exactly. For every pure input rho^tensor4 the output is |Theta psi tensor psi^tensor3><Theta psi tensor psi^tensor3|; it is positive, with trace1. Convexity gives positivity for rotationally averaged product moments.',scope='Checks the complete 1|3 transpose/time-reversal formula and half-integer signs. This validates the positive-map constraint; it does not itself prove nonexistence.')
(HERE/'positive-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
