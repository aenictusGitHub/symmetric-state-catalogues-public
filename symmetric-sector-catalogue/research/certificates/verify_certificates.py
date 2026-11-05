"""Verify every supplied rational certificate, without a numerical optimizer.

Run: python3 verify_certificates.py
Only Python's standard library is required.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from spectral_bounds import purity_eigenvalue, sector_eigenvalue
from certificate_basis import constraints, multipole_eigenvalue

rows=json.loads(Path(__file__).with_name('lp-bounds.json').read_text())
for row in rows:
 n,k,mode=row['n'],row['k'],row['mode'];B=F(row['value'])
 basis=dict(constraints(n));ds=range(0,n+1,2)
 coeffs={name:F(value) for name,value in row['certificate'].items()}
 assert all(value>=0 for value in coeffs.values()),row
 for idx,d in enumerate(ds):
  s=sector_eigenvalue(n,k,d)
  target=B-s if mode=='max' else s-B
  obtained=sum((value*basis[name][idx] for name,value in coeffs.items()),F(0))
  assert obtained==target,(n,k,mode,n-d,obtained,target)
 # The relaxation extremal vertex must satisfy all exact constraints too.
 p=[F(row['pJ_vertex'][str(n-d)]) for d in ds]
 assert sum(p)==1
 for name,g in basis.items():assert sum(x*y for x,y in zip(p,g))>=0,(row,name)
 assert sum(p[i]*sector_eigenvalue(n,k,d) for i,d in enumerate(ds))==B
# Cross-check universal relations for the spectral basis.
for n in range(2,13):
 for d in range(0,n+1,2):
  assert purity_eigenvalue(n,d,0)==purity_eigenvalue(n,d,n)==1
  assert all(purity_eigenvalue(n,d,t)==purity_eigenvalue(n,d,n-t) for t in range(n+1))
  assert sum(sector_eigenvalue(n,k,d) for k in range(n+1))==2**n
  assert sum(multipole_eigenvalue(n,L,d) for L in range(n+1))==1
print(f'PASS: {len(rows)} exact rational bounds, all nonnegative coefficients, all identities, all relaxation vertices.')
print('This proves bounds. Attainment requires a separately verified physical state.')
