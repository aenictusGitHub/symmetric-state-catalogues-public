"""Exploratory unrestricted AC2 search. Numerical residuals are not proofs."""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent

def problem(n,k):
 d=n+1;j=n/2
 z=np.diag(j-np.arange(d))/j
 plus=np.diag(np.sqrt(np.arange(1,d)*(d-np.arange(1,d))),1)/j
 x=(plus+plus.T)/2;y=(plus-plus.T)/(2j)
 c=(j+1)/(3*j)
 ops=np.array([np.eye(d),x,y,z,x@x-y@y,x@y+y@x,x@z+z@x,y@z+z@y,z@z-c*np.eye(d)])
 ii=np.arange(k);aa,bb=np.triu_indices(k,1);eye=np.eye(k)
 def pack(A):return np.concatenate([A[:,ii,ii].real,np.sqrt(2)*A[:,aa,bb].real,np.sqrt(2)*A[:,aa,bb].imag],axis=1).reshape(-1)
 def unpack(v):return v[:d*k].reshape(d,k)+1j*v[d*k:].reshape(d,k)
 def fun(v):
  B=unpack(v);C=np.einsum('ra,lrs,sb->lab',B.conj(),ops,B);C[0]-=eye;return pack(C)
 def jac(v):
  B=unpack(v);Z=ops@B
  row=np.einsum('ac,lrb->labrc',eye,Z);col=np.einsum('bc,lra->labrc',eye,Z.conj())
  deriv=np.concatenate([(row+col).reshape(9,k,k,-1),(1j*(col-row)).reshape(9,k,k,-1)],axis=-1)
  return np.concatenate([deriv[:,ii,ii].real,np.sqrt(2)*deriv[:,aa,bb].real,np.sqrt(2)*deriv[:,aa,bb].imag],axis=1).reshape(9*k*k,-1)
 return fun,jac,unpack

def run(n,k,restarts=8,seed=581,steps=1000):
 fun,jac,unpack=problem(n,k);rng=np.random.default_rng(seed);best=1e9;start=time.time()
 for r in range(restarts):
  B=np.linalg.qr(rng.normal(size=(n+1,k))+1j*rng.normal(size=(n+1,k)))[0]
  v=np.r_[B.real.reshape(-1),B.imag.reshape(-1)]
  if r==0:
   direction=rng.normal(size=len(v));eps=1e-6
   assert np.max(abs((fun(v+eps*direction)-fun(v-eps*direction))/(2*eps)-jac(v)@direction))<1e-7
  sol=least_squares(fun,v,jac=jac,max_nfev=steps,ftol=1e-12,gtol=1e-12,xtol=1e-12)
  val=float(np.linalg.norm(sol.fun));print(n,k,r,sol.nfev,val,round(time.time()-start,1),flush=True)
  if val<best:
   best=val;B=unpack(sol.x)
   (HERE/f'general-n{n}-k{k}.json').write_text(json.dumps({'n':n,'k':k,'residual':val,'attempt':r,'seed':seed,'basis':[[[float(z.real),float(z.imag)] for z in B[:,i]] for i in range(k)],'note':'Exploratory numerical search only; neither existence nor nonexistence is certified.'},indent=2)+'\n')
  if best<1e-10:break
 return best
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=8);p.add_argument('--steps',type=int,default=1000);a=p.parse_args();run(a.n,a.k,restarts=a.restarts,steps=a.steps)
