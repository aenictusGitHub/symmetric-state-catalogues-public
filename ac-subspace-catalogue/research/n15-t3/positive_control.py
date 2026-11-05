#!/usr/bin/env python3
"""Supplementary N13 k2 AC3 positive control for the Young22 moment equations.
Requires numpy, scipy, and sympy. Its numerical rotation is a diagnostic,
not a premise of the N15 exact exclusion certificate.
"""
import json,numpy as np,sympy as s
from scipy.linalg import expm
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent
D=json.loads((HERE.parents[1]/'data/catalogue.json').read_text())
sub=next(x for x in D['subspaces']if x['id']=='radical-n13-k2-t3');n=sub['n'];d=n+1
U=s.Matrix.hstack(*[s.Matrix([s.sympify(x)for x in b['exactCoefficients']])for b in sub['basis']])
assert U.shape==(14,2)
# Full-span cross conditions, not only separate anticoherence of each vector.
z_exact=s.diag(*[s.Rational(2*i-n,2)for i in range(d)])
up_exact=s.zeros(d)
for i in range(d-1):up_exact[i+1,i]=s.sqrt((n-i)*(i+1))
compression_checks=[]
for q in range(4):
 for p in range(4-q):
  matrix=U.conjugate().T*(up_exact**q)*(z_exact**p)*U
  target=(sum(z_exact[i,i]**p for i in range(d))/d)*s.eye(2) if q==0 else s.zeros(2)
  assert all(s.simplify(x)==0 for x in matrix-target)
  compression_checks.append({'raisingPower':q,'diagonalPower':p,'exactCrossCompression':True})
# First-leg contractions of the normalized Slater determinant, evaluated directly.
W=(U[:,0]*U[:,1].T-U[:,1]*U[:,0].T)/s.sqrt(2)
rows=[]
for ell in range(1,4):
 f=[(2*i-n)if ell==1 else 3*(2*i-n)**2-n*(n+2)if ell==2 else(2*i-n)*(5*(2*i-n)**2-3*n*(n+2)+4)for i in range(d)]
 C=W.T*s.diag(*f)*W.conjugate();nz=[s.simplify(x)for x in C if x!=0];assert all(x==0 for x in nz)
 rows.append({'ell':ell,'allFirstLegContractionEntriesExactlyZero':True})
# Rotate the exact complex subspace and independently verify physical tensor coordinates.
V=np.array(U.evalf(18),complex);z=np.diag(np.arange(d)-n/2);p=np.zeros((d,d))
for i in range(d-1):p[i+1,i]=np.sqrt((n-i)*(i+1))
x=(p+p.T)/2;y=(p-p.T)/2j
V=expm(-1j*(.37*x+.91*y+.23*z))@V
w=(np.outer(V[:,0],V[:,1])-np.outer(V[:,1],V[:,0]))/np.sqrt(2)
plucker=max(abs(w[a,b]*w[c,e]-w[a,c]*w[b,e]+w[a,e]*w[b,c])for a,b,c,e in combinations(range(d),4))
errors=[]
for ell in range(1,4):
 f=np.array([(2*i-n)if ell==1 else 3*(2*i-n)**2-n*(n+2)if ell==2 else(2*i-n)*(5*(2*i-n)**2-3*n*(n+2)+4)for i in range(d)],float)
 errors.append(float(np.linalg.norm(w.T@np.diag(f)@w.conj())/max(abs(f))))
assert max(errors)<1e-14 and plucker<1e-14
report={'passed':True,'sourceSubspaceId':sub['id'],'n':13,'t':3,'k':2,'exactCrossCompressionChecks':compression_checks,'exactFirstLegContractions':rows,'rotatedComplexBasisNormError':float(np.linalg.norm(V.conj().T@V-np.eye(2))),'rotatedWedgeSquaredNorm':float(np.linalg.norm(w)**4),'maximumRotatedPluckerResidual':float(plucker),'rotatedFirstLegNormalizedResiduals':errors,'scope':'Positive control for the necessary Young22 support and first-leg equations; not an N15 nonexistence proof.'}
(HERE/'positive-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
