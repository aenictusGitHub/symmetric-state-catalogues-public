#!/usr/bin/env python3
"""Bounded complex Grassmann search; numerical evidence only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
N=13; d=14;k=3;j=N/2
m=j-np.arange(d)
Jz=np.diag(m).astype(complex)
Jp=np.zeros((d,d),complex)
for r in range(1,d): Jp[r-1,r]=np.sqrt(r*(N-r+1))
Jx=(Jp+Jp.conj().T)/2; Jy=(Jp-Jp.conj().T)/(2j)
raw=[Jx,Jy,Jz,Jx@Jx-Jy@Jy,(2*Jz@Jz-Jx@Jx-Jy@Jy)/np.sqrt(3),Jx@Jy+Jy@Jx,Jx@Jz+Jz@Jx,Jy@Jz+Jz@Jy]
T=np.array([np.eye(d)]+[a/np.linalg.norm(a) for a in raw])
target=np.zeros((9,k,k),complex);target[0]=np.eye(k)
inds=np.triu_indices(k,1)
sqrt2=np.sqrt(2)
def pack(C):
 return np.concatenate([C.diagonal(axis1=-2,axis2=-1).real,C[...,inds[0],inds[1]].real*sqrt2,C[...,inds[0],inds[1]].imag*sqrt2],axis=-1)
def unpack(x): return (x[:d*k]+1j*x[d*k:]).reshape(d,k)
def fun(x):
 B=unpack(x); return pack(np.einsum('ri,arj->aij',B.conj(),T@B)-target).ravel()
def jac(x):
 B=unpack(x); AB=T@B
 # Complex derivative arrays [operator, r, a, i, j].
 eye=np.eye(k)
 first=np.einsum('ai,orj->oraij',eye,AB)
 second=np.einsum('aj,ori->oraij',eye,AB.conj())
 dx=pack(first+second).transpose(0,3,1,2).reshape(81,42)
 dy=pack(-1j*first+1j*second).transpose(0,3,1,2).reshape(81,42)
 return np.concatenate((dx,dy),axis=1)
def coords(B):return np.r_[B.real.ravel(),B.imag.ravel()]
rng=np.random.default_rng(130203)
B0=np.linalg.qr(rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)))[0]
x0=coords(B0); direction=rng.normal(size=84); direction/=np.linalg.norm(direction)
eps=1e-6
jac_error=float(np.max(np.abs((fun(x0+eps*direction)-fun(x0-eps*direction))/(2*eps)-jac(x0)@direction)))
assert jac_error<1e-8,jac_error
seeds=[('random-0',B0)]
p=Path(__file__).resolve().parent.parent/'algebraic-k3/n14/basis.json'
source=json.loads(p.read_text())['subspaces'][0]
B14=np.array([[complex(*c) for c in row['coefficients']] for row in source['basis']]).T
z14=np.diag(7-np.arange(15)); plus14=np.zeros((15,15),complex)
for r in range(1,15): plus14[r-1,r]=np.sqrt(r*(15-r))
x14=(plus14+plus14.conj().T)/2;y14=(plus14-plus14.conj().T)/(2j)
raw14=[x14,y14,z14,x14@x14-y14@y14,(2*z14@z14-x14@x14-y14@y14)/np.sqrt(3),x14@y14+y14@x14,x14@z14+z14@x14,y14@z14+z14@y14]
control=max(float(np.linalg.norm(B14.conj().T@(a/np.linalg.norm(a))@B14)) for a in raw14)
assert control<1e-12,control
# The two one-qubit contractions of the exact N14 code provide nearby physical seeds.
for name,B in [('n14-zero-contraction',np.sqrt((14-np.arange(14))/14)[:,None]*B14[:14]),('n14-one-contraction',np.sqrt((np.arange(14)+1)/14)[:,None]*B14[1:])]:
 seeds.append((name,np.linalg.qr(B)[0]))
for a in range(1,7): seeds.append((f'random-{a}',np.linalg.qr(rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)))[0]))
t0=time.monotonic();deadline=t0+110; results=[];best=(float('inf'),None,None)
class TimeLimit(Exception):pass
for name,B in seeds:
 if time.monotonic()>deadline:break
 localbest=[np.linalg.norm(fun(coords(B))),coords(B).copy()]
 def fbounded(x):
  if time.monotonic()>deadline: raise TimeLimit()
  r=fun(x); z=np.linalg.norm(r)
  if z<localbest[0]: localbest[:]=[z,x.copy()]
  return r
 try:
  opt=least_squares(fbounded,coords(B),jac=jac,method='trf',max_nfev=220,ftol=1e-11,xtol=1e-11,gtol=1e-11)
  status=int(opt.status);nfev=int(opt.nfev)
 except TimeLimit:status='global-time-limit';nfev=None
 norm,x=localbest;B=unpack(x)
 item={'seed':name,'residualL2':float(norm),'residualMax':float(np.max(np.abs(fun(x)))),'status':status,'nfev':nfev}
 results.append(item)
 if norm<best[0]: best=(norm,x.copy(),name)
 print(json.dumps(item),flush=True)
 if norm<1e-10:break
B=unpack(best[1]);gram=B.conj().T@B
# A QR step tests the physical compression independently of the residual parameterization.
U=np.linalg.qr(B)[0]
physical=[float(np.linalg.norm(U.conj().T@a@U)) for a in T[1:]]
report={'n':13,'dimension':3,'order':2,'status':'numerical-candidate-only' if max(physical)<1e-9 else 'no-verified-three-dimensional-construction-found','seedCount':len(results),'seedGenerator':130203,'elapsedSeconds':time.monotonic()-t0,'analyticJacobianDirectionalError':jac_error,'constraints':81,'realVariables':84,'tensorNormalization':'Hilbert-Schmidt unit norm for eight rank-one and rank-two Hermitian tensors; identity unscaled','restarts':results,'bestSeed':best[2],'bestResidualL2':float(best[0]),'bestGramResidual':float(np.linalg.norm(gram-np.eye(3))),'orthonormalizedTensorCompressionNorms':physical,'bestBasisAfterQR':[[[float(z.real),float(z.imag)] for z in U[:,a]] for a in range(3)],'interpretation':'A failed finite numerical search does not prove nonexistence. A small residual alone would not be an exact existence certificate.'}
report['positiveControl']={'n':14,'dimension':3,'knownExactCodeMaxCompressionNorm':control}
(Path(__file__).resolve().parent/'direct-search.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({key:report[key] for key in ['status','seedCount','elapsedSeconds','bestSeed','bestResidualL2','orthonormalizedTensorCompressionNorms']}),flush=True)
