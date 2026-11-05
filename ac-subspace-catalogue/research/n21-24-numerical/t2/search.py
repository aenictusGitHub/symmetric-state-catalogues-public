#!/usr/bin/env python3
"""Reproducible full-complex AC2 multistart searches; no nonexistence inference."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('MKL_NUM_THREADS','1')
from pathlib import Path
import argparse,json,time,math,hashlib
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent;CAT=HERE.parents[2]

def problem(n,k,indices=None):
 d=n+1;j=n/2;indices=np.arange(d)if indices is None else np.array(indices);d0=len(indices)
 z=np.diag(j-np.arange(d))/j;p=np.diag(np.sqrt(np.arange(1,d)*(d-np.arange(1,d))),1)/j
 x=(p+p.T)/2;y=(p-p.T)/(2j);c=(j+1)/(3*j)
 ops=np.array([np.eye(d),x,y,z,x@x-y@y,x@y+y@x,x@z+z@x,y@z+z@y,z@z-c*np.eye(d)])[:,indices][:,:,indices]
 ii=np.arange(k);aa,bb=np.triu_indices(k,1);eye=np.eye(k)
 def pack(A):return np.concatenate([A[:,ii,ii].real,np.sqrt(2)*A[:,aa,bb].real,np.sqrt(2)*A[:,aa,bb].imag],axis=1).reshape(-1)
 def unpack(v):return v[:d0*k].reshape(d0,k)+1j*v[d0*k:].reshape(d0,k)
 def fun(v):
  B=unpack(v);C=np.einsum('ra,lrs,sb->lab',B.conj(),ops,B);C[0]-=eye;return pack(C)
 def jac(v):
  B=unpack(v);Z=ops@B;row=np.einsum('ac,lrb->labrc',eye,Z);col=np.einsum('bc,lra->labrc',eye,Z.conj())
  deriv=np.concatenate([(row+col).reshape(9,k,k,-1),(1j*(col-row)).reshape(9,k,k,-1)],axis=-1)
  return np.concatenate([deriv[:,ii,ii].real,np.sqrt(2)*deriv[:,aa,bb].real,np.sqrt(2)*deriv[:,aa,bb].imag],axis=1).reshape(9*k*k,-1)
 return fun,jac,unpack,indices

def coords(B):return np.r_[B.real.ravel(),B.imag.ravel()]
def basis_json(B):return [{'coefficients':[[float(z.real),float(z.imag)]for z in B[:,i]]}for i in range(B.shape[1])]
def load_basis(path):
 d=json.loads(path.read_text());return np.array([[complex(*z)for z in row['coefficients']]for row in d['basis']]).T

def check(B,n,t=2):
 mats=[np.array([[B[a+h,i]*math.sqrt(math.comb(t,a)*math.comb(n-t,h)/math.comb(n,a+h))for h in range(n-t+1)]for a in range(t+1)])for i in range(B.shape[1])]
 errs=[float(np.linalg.norm(A@D.conj().T-(np.eye(t+1)/(t+1)if a==b else np.zeros((t+1,t+1)))))for a,A in enumerate(mats)for b,D in enumerate(mats)]
 return {'orthonormalityFrobenius':float(np.linalg.norm(B.conj().T@B-np.eye(B.shape[1]))),'maxOrderedCrossReductionFrobenius':max(errs),'orderedReducedMatrixCount':len(errs),'complexReducedEntriesChecked':len(errs)*(t+1)**2,'allOrderedCrossReductionErrors':errs}

def sources(n,k):
 result=[]
 if k==4:
  p=CAT/'research/n20-t2/basis.json';B=load_basis(p)
  for shift in range(n-20+1):
   A=np.zeros((n+1,k),complex);A[shift:shift+21]=B;result.append((f'exact-n20-embedded-shift-{shift}',np.linalg.qr(A)[0]))
 if k==5:
  p=CAT/'research/t2-general/n24-k5/basis.json';B=load_basis(p)
  for a in range(24-n+1):
   A=np.array([B[r+a]*math.sqrt(math.comb(24-n,a)*math.comb(n,r)/math.comb(24,r+a))for r in range(n+1)]);result.append((f'exact-n24-{24-n}-qubit-contraction-{a}',np.linalg.qr(A)[0]))
 return result

def run(n,k,restarts=12,steps=1000,seed=20261103,parity=None):
 indices=np.arange(n+1)if parity is None else np.arange(parity,n+1,2);fun,jac,unpack,indices=problem(n,k,indices);rng=np.random.default_rng(seed+n*101+k*7);d=len(indices)
 seeds=[]
 if parity is None:seeds=sources(n,k)
 seeds+=[(f'complex-gaussian-{i}',np.linalg.qr(rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)))[0])for i in range(restarts)]
 best=np.inf;attempts=[];start=time.monotonic();bestB=None;bestname=None
 stem=f'n{n}-k{k}'+(''if parity is None else f'-parity{parity}');out=HERE/(stem+'.json')
 for index,(name,B)in enumerate(seeds):
  v=coords(B)
  if index==0:
   q=rng.normal(size=len(v));q/=np.linalg.norm(q);h=1e-6;derivative_error=float(np.max(abs((fun(v+h*q)-fun(v-h*q))/(2*h)-jac(v)@q)));assert derivative_error<1e-8
  before=time.monotonic();sol=least_squares(fun,v,jac=jac,max_nfev=steps,ftol=3e-13,gtol=3e-13,xtol=3e-13)
  C=np.zeros((n+1,k),complex);C[indices]=unpack(sol.x);C=np.linalg.qr(C)[0];physical=check(C,n)
  residual=float(np.linalg.norm(sol.fun));item={'seedName':name,'nfev':int(sol.nfev),'optimizerStatus':int(sol.status),'constraintResidualL2':residual,'maxCrossReductionFrobenius':physical['maxOrderedCrossReductionFrobenius'],'orthonormalityFrobenius':physical['orthonormalityFrobenius'],'elapsedSeconds':time.monotonic()-before}
  attempts.append(item)
  if physical['maxOrderedCrossReductionFrobenius']<best:best=physical['maxOrderedCrossReductionFrobenius'];bestB=C;bestname=name;bestcheck=physical
  report={'n':n,'t':2,'dimension':k,'id':f'numerical-n{n}-k{k}-t2','evidence':'numerical','found':best<1e-10,'searchRestriction':'arbitrary complex Dicke coefficients'if parity is None else f'single Dicke parity {parity}','randomSeed':seed+n*101+k*7,'analyticJacobianDirectionalError':derivative_error,'attempts':attempts,'bestSeed':bestname,'bestCrossReductionFrobenius':best,'basis':basis_json(bestB),'independentPhysicalCheck':bestcheck,'totalElapsedSeconds':time.monotonic()-start,'interpretation':'A successful floating-point witness is numerical evidence only. Failure of finite searches is not a nonexistence proof. Dimension maximality is not inferred from search failures.'}
  out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'n':n,'k':k,'parity':parity,**item}),flush=True)
  if best<1e-12:break
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=12);p.add_argument('--steps',type=int,default=1000);p.add_argument('--seed',type=int,default=20261103);p.add_argument('--parity',type=int);a=p.parse_args();run(a.n,a.k,a.restarts,a.steps,a.seed,a.parity)
