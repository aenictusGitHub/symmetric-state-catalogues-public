#!/usr/bin/env python3
"""Numerical convention controls for the exact spin-three inertia witness.
These floating-point checks are diagnostics, not premises of the proof.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def generators(j):
 d=int(2*j+1);p=np.zeros((d,d),complex)
 for r in range(1,d):p[r-1,r]=np.sqrt(r*(2*j-r+1))
 return [(p+p.conj().T)/2,(p-p.conj().T)/(2j),np.diag(j-np.arange(d))]
def verify():
 raw=(ROOT/'data/catalogue.json').read_bytes();data=json.loads(raw);J=generators(9);S=generators(3);X=sum(np.kron(a,b)for a,b in zip(J,S));I=np.eye(133)
 W=X@X@X@X-2*X@X@X-653*X@X-1146*X+8640*I;scale=319320
 values=np.linalg.eigvalsh(W);expected=np.concatenate([np.full(2*L+1,v)for L,v in zip(range(6,13),[319320,-6264,-63720,-6264,-6264,-111600,-6264])]);error=float(np.max(abs(np.sort(values)-np.sort(expected))))
 assert error<1e-8 and np.sum(values>0)==13 and np.sum(values<0)==120
 partial=np.einsum('aiaj->ij',W.reshape(19,7,19,7));assert np.linalg.norm(partial)<1e-8
 rotations=[np.eye(19),expm(-.381j*J[2])@expm(.917j*J[1])@expm(-.233j*J[0])]
 records=[]
 for sid,t in [('line-acliege-n18-t5',5),('family51-n18',4),('radical-n18-k2-t3',3)]:
  source=next(s for s in data['subspaces']if s['id']==sid);B=np.array([[complex(*z)for z in v['coefficients']]for v in source['basis']]).T
  norms=[]
  for rot in rotations:
   F=np.kron(rot@B,np.eye(7));norms.append(float(np.linalg.norm(F.conj().T@W@F)/scale))
  if t>=4:assert max(norms)<2e-14
  else:assert min(norms)>1e-5
  records.append({'subspaceId':sid,'verifiedACOrder':t,'dimension':B.shape[1],'normalizedCompressionNormsOriginalAndRotated':norms,'compressionExpectedZero':t>=4})
 return {'passed':True,'n':18,'order':4,'auxiliarySpin':3,'polynomial':'X^4 - 2X^3 - 653X^2 - 1146X + 8640I','physicalDimension':133,'positiveInertia':13,'negativeInertia':120,'maximumEigenvalueError':error,'partialTraceNorm':float(np.linalg.norm(partial)),'controls':records,'scope':'Numerical spin/normalization controls only. Exact upper and lower proofs are verified separately.'}
if __name__=='__main__':
 r=verify();(HERE/'positive-control.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
