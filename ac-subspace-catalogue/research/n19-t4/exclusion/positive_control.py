#!/usr/bin/env python3
"""Exact AC4-plane positive control for the Young22 physical contraction equations.
The N=60 control is a separate analytical construction. Exact sparse moment
checks prove AC4 of its entire two-dimensional span. A generic complex rotation
then checks normalization, Pluecker support and direct first-leg contractions.
Numerical checks are diagnostics only and are not exclusion-certificate premises.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, lcm, sqrt
from itertools import combinations
import json
import numpy as np
from scipy.linalg import expm


def verify():
 n=60; d=n+1; j=F(n,2)
 # These are probabilities of each signed magnetic point, not the pair total.
 distributions=[{0:F(6667,50625),-15:F(53041,151875),15:F(53041,151875),
                 -30:F(12896,151875),30:F(12896,151875)},
                {-5:F(889,7500),5:F(889,7500),-10:F(617,3750),10:F(617,3750),
                 -25:F(1627,7500),25:F(1627,7500)}]
 expected=[sum((F(i)-j)**p for i in range(d))/d for p in range(5)]
 assert expected==[1,0,310,0,172918]
 support=sorted(m for dist in distributions for m in dist)
 assert len(support)==len(set(support)) and all(b-a>=5 for a,b in zip(support,support[1:]))
 for dist in distributions:
  assert all(w>0 for w in dist.values())
  assert [sum(w*F(m)**p for m,w in dist.items())for p in range(5)]==expected
 # Disjoint supports kill zero-shift cross moments. Union gaps5 kill every
 # nonzero magnetic shift through rank4, including cross-vector matrix elements.
 U=np.zeros((d,2),complex)
 for col,dist in enumerate(distributions):
  for m,w in dist.items(): U[int(m+j),col]=sqrt(float(w))
 z=np.diag(np.arange(d)-float(j)); up=np.zeros((d,d))
 for i in range(n): up[i+1,i]=sqrt((n-i)*(i+1))
 x=(up+up.T)/2; y=(up-up.T)/2j
 V=expm(-1j*(.37*x+.91*y+.23*z))@U
 # Unit-length wedge gives unit-length w tensor w in the Young(2,2) summand.
 w=(np.outer(V[:,0],V[:,1])-np.outer(V[:,1],V[:,0]))/sqrt(2)
 plucker=max(abs(w[a,b]*w[c,e]-w[a,c]*w[b,e]+w[a,e]*w[b,c])
             for a,b,c,e in combinations(range(d),4))
 assert plucker<2e-14
 polys=[]; errors=[]
 for ell in range(5):
  f=[F((2*i-n)**ell)for i in range(d)]
  for p in polys:
   q=sum(a*b for a,b in zip(f,p))/sum(a*a for a in p)
   f=[a-q*b for a,b in zip(f,p)]
  den=lcm(*(a.denominator for a in f)); f=[int(a*den)for a in f]; common=gcd(*f)
  f=[a//common for a in f]; polys.append(f)
  assert all(((n-i)*(i+1)*(f[i]-f[i+1]) if i<n else 0)
             +(i*(n-i+1)*(f[i]-f[i-1]) if i>0 else 0)==ell*(ell+1)*f[i]for i in range(d))
  if ell:
   fn=np.array(f,float)/max(abs(a)for a in f)
   contraction=w.T@np.diag(fn)@w.conj()
   err=float(np.linalg.norm(contraction)); assert err<2e-14
   errors.append({'rank':ell,'normalizedFirstLegContractionNorm':err,
                  'rankFourTensorContractionFrobeniusNorm':err*float(np.linalg.norm(w)**2)})
 normerr=float(np.linalg.norm(V.conj().T@V-np.eye(2)))
 wedge_norm=float(np.linalg.norm(w)**2); moment_trace=wedge_norm**2
 assert normerr<2e-14 and abs(moment_trace-1)<2e-14
 return {'passed':True,'n':n,'t':4,'dimension':2,
         'exactSignedMagneticPointProbabilities':[{str(k):str(v)for k,v in dist.items()}for dist in distributions],
         'exactDiagonalMomentsThroughFour':[str(x)for x in expected],
         'allProbabilitiesStrictlyPositive':True,'fullSpanAC4CrossConditionsVerifiedExactly':True,
         'exactBasisGramMatrix':'identity(2)','exactNonzeroNormalizedYoung22MomentTrace':'1',
         'disjointSupports':True,'minimumUnionSupportGap':min(b-a for a,b in zip(support,support[1:])),
         'rotatedComplexBasisGramError':normerr,'normalizedWedgeSquaredNorm':wedge_norm,
         'normalizedNonzeroYoung22MomentTrace':moment_trace,
         'maximumRotatedPluckerResidual':float(plucker),'rotatedFirstLegContractions':errors,
         'scope':'A genuine normalized nonzero AC4-plane Young22 moment satisfies the necessary support and first-leg equations. The rotation is a supplementary numerical diagnostic, not a premise of either exact exclusion proof.'}

if __name__=='__main__':
 report=verify(); root=Path(__file__).resolve().parent
 (root/'positive-control.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
