#!/usr/bin/env python3
"""Reproducible multistart numerical search; results are not optimality proofs."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, check_grad

ROOT=Path(__file__).resolve().parent/'runs'
ROOT.mkdir(exist_ok=True)
CAT=Path(__file__).resolve().parents[1]/'catalogue-before-update.json'

class Sector:
 def __init__(self,n,k):
  self.n,self.k=n,k
  self.w=np.zeros(n//2+1)
  for r in range(k+1):self.w[min(r,n-r)]+=math.comb(n,k)*(-1)**(k-r)*math.comb(k,r)*2**r
  self.parts=[]
  for r,w in enumerate(self.w):
   if not w:continue
   aa,bb=np.indices((r+1,n-r+1)); ids=aa+bb
   fac=np.array([[math.sqrt(math.comb(r,a)*math.comb(n-r,b)/math.comb(n,a+b)) for b in range(n-r+1)] for a in range(r+1)])
   self.parts.append((w,ids,fac))
 def value_grad(self,x):
  norm=np.linalg.norm(x); y=x/norm
  c=y[:self.n+1]+1j*y[self.n+1:]
  value=0.; g=np.zeros(self.n+1,complex)
  for w,ids,fac in self.parts:
   m=c[ids]*fac; rho=m@m.conj().T
   value+=w*np.sum(np.abs(rho)**2)
   a=(rho@m)*fac
   g+=4*w*(np.bincount(ids.ravel(),weights=a.real.ravel(),minlength=self.n+1)+1j*np.bincount(ids.ravel(),weights=a.imag.ravel(),minlength=self.n+1))
  grad=np.r_[g.real,g.imag];grad=(grad-y*np.dot(grad,y))/norm
  return float(value),grad
 def sectors(self,c):
  x=np.r_[c.real,c.imag]
  return [Sector(self.n,k).value_grad(x)[0] for k in range(self.n+1)]

def vectors(n,rng):
 out=[]
 for j in range(n+1):
  c=np.zeros(n+1,complex);c[j]=1;out.append(c)
 c=np.zeros(n+1,complex);c[[0,-1]]=1/math.sqrt(2);out.append(c)
 for e in json.load(CAT.open())['entries']:
  if e['n']==n and e.get('coefficients'):
   try:
    c=np.array([complex(*a) if isinstance(a,list) else complex(a) for a in e['coefficients']]);c/=np.linalg.norm(c);out.append(c)
   except:pass
 return out

def search(n,k,mode,starts=32,seed=20261028):
 rng=np.random.default_rng(seed+n*100+k*2+(mode=='max'))
 obj=Sector(n,k);sgn=-1 if mode=='max' else 1
 seeds=vectors(n,rng)
 # Pure Dicke vectors are stationary: perturb them so optimizer can escape.
 seeds=[c+0.03*(rng.normal(size=n+1)+1j*rng.normal(size=n+1)) for c in seeds]
 seeds+= [rng.normal(size=n+1)+1j*rng.normal(size=n+1) for _ in range(starts)]
 best=None; history=[];start=time.time()
 def f(x):
  value,grad=obj.value_grad(x);return sgn*value,sgn*grad
 for j,c in enumerate(seeds):
  c/=np.linalg.norm(c);x=np.r_[c.real,c.imag]
  res=minimize(f,x,jac=True,method='L-BFGS-B',options={'maxiter':1200,'ftol':1e-14,'gtol':2e-8,'maxls':40})
  value=sgn*res.fun
  history.append({'start':j,'value':float(value),'success':bool(res.success),'nit':int(res.nit),'gradientNorm':float(np.linalg.norm(res.jac))})
  if best is None or sgn*value<sgn*best['value']:
   y=res.x/np.linalg.norm(res.x);c=y[:n+1]+1j*y[n+1:];c*=np.exp(-1j*np.angle(c[np.argmax(np.abs(c))]))
   best={'n':n,'k':k,'mode':mode,'value':float(value),'coefficients':[[float(a.real),float(a.imag)] for a in c], 'sectorVector':obj.sectors(c),'gradientNorm':float(np.linalg.norm(res.jac)),'optimizerSuccess':bool(res.success),'optimizerMessage':str(res.message)}
 best.update(seed=seed,starts=len(seeds),history=history,elapsed=time.time()-start,status='numerical',warning='Numerical candidate; not a proof of global optimality.')
 (ROOT/f'n{n}-s{k}-{mode}.json').write_text(json.dumps(best,indent=2)+'\n')
 print(f'n{n} S{k} {mode}: {best["value"]:.12g}; starts={len(seeds)}; {best["elapsed"]:.1f}s',flush=True)
 return best

if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='check':
  obj=Sector(7,5);x=np.random.default_rng(4).normal(size=16)
  print('gradient error',check_grad(lambda x:obj.value_grad(x)[0],lambda x:obj.value_grad(x)[1],x))
 elif len(sys.argv)>1:
  search(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],int(sys.argv[4]) if len(sys.argv)>4 else 32)
 else:
  tasks=[(10,8,'max'),(10,9,'max')]+[(11,k,'max') for k in range(6,11)]+[(12,k,'max') for k in range(6,12)]+[(9,6,'min')]+[(n,4,'min') for n in (10,11,12)]
  results=[search(*task) for task in tasks]
  (ROOT/'all-results.json').write_text(json.dumps(results,indent=2)+'\n')
