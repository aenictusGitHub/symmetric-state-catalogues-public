#!/usr/bin/env python3
"""Reproducible full-complex multistart pure-AC search; numerical evidence only."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ.setdefault(key,'1')
from pathlib import Path
import argparse,json,time,math,hashlib
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
SOURCES={21:4,22:5,23:5,24:7}
def load(n):
 a=np.loadtxt(HERE/f'inputs/acsdk.{n}.{SOURCES[n]}.txt');v=a[:,0]+1j*a[:,1];assert len(v)==n+1
 return v/np.linalg.norm(v)
def ops(n,t):
 d=n+1;Jp=np.diag(np.sqrt(np.arange(1,d)*(d-np.arange(1,d))),1).astype(complex);Jm=Jp.conj().T
 out=[np.eye(d)/np.sqrt(d)]
 for ell in range(1,t+1):
  T=np.linalg.matrix_power(Jp,ell);T/=np.linalg.norm(T);qs={ell:T}
  for q in range(ell,0,-1):qs[q-1]=(Jm@qs[q]-qs[q]@Jm)/math.sqrt((ell+q)*(ell-q+1))
  out.append(qs[0])
  for q in range(1,ell+1):out.extend([(qs[q]+qs[q].conj().T)/np.sqrt(2),(qs[q]-qs[q].conj().T)/(1j*np.sqrt(2))])
 out=np.array(out);gram=np.einsum('aij,bij->ab',out.conj(),out).real
 assert np.max(abs(gram-np.eye(len(out))))<1e-9
 assert np.max(abs(out-out.conj().transpose(0,2,1)))<1e-10
 return out

def reductions(v,n,t):
 rows=[]
 for order in range(1,t+1):
  M=np.array([[v[a+b]*math.sqrt(math.comb(order,a)*math.comb(n-order,b)/math.comb(n,a+b))for b in range(n-order+1)]for a in range(order+1)])
  E=M@M.conj().T-np.eye(order+1)/(order+1)
  rows.append({'order':order,'frobeniusResidual':float(np.linalg.norm(E)),'maxEntryResidual':float(np.max(abs(E)))})
 return rows

def search(n,t,restarts,steps,seed):
 d=n+1;T=ops(n,t);target=np.zeros(len(T));target[0]=1/np.sqrt(d);rng=np.random.default_rng(seed)
 def pack(v):return np.r_[v.real,v.imag]
 def unpack(x):return x[:d]+1j*x[d:]
 def fun(x):
  v=unpack(x);return np.einsum('i,ai->a',v.conj(),T@v).real-target
 def jac(x):
  y=T@unpack(x);return 2*np.c_[y.real,y.imag]
 v=load(n);x=pack(v);dx=rng.normal(size=2*d);eps=1e-6
 jacerr=float(np.max(abs((fun(x+eps*dx)-fun(x-eps*dx))/(2*eps)-jac(x)@dx)));assert jacerr<1e-8
 attempts=[];best=None;beg=time.monotonic()
 for r in range(restarts):
  if r==0:v=load(n);kind='local-library'
  elif r<4:
   v=load(n)+[.01,.1,.5][r-1]*(rng.normal(size=d)+1j*rng.normal(size=d));kind='perturbed-library'
  elif r<8:
   high=load(24);offset=(r-4)%(25-d+1);v=high[offset:offset+d];kind='truncated-n24-ac7'
  else:v=rng.normal(size=d)+1j*rng.normal(size=d);kind='random-complex'
  v/=np.linalg.norm(v);tic=time.monotonic()
  sol=least_squares(fun,pack(v),jac=jac,method='trf',max_nfev=steps,ftol=2e-14,gtol=2e-14,xtol=2e-14)
  v=unpack(sol.x);rawnorm=float(np.linalg.norm(v));v/=rawnorm
  # Fix one global phase for export only.
  phase=v[np.argmax(abs(v))];v*=phase.conjugate()/abs(phase)
  checks=reductions(v,n,t);err=max(x['frobeniusResidual']for x in checks)
  row={'attempt':r,'seed':seed,'initialization':kind,'functionEvaluations':int(sol.nfev),'jacobianEvaluations':int(sol.njev),'solverStatus':int(sol.status),'tensorResidualL2':float(np.linalg.norm(sol.fun)),'maximumReductionFrobeniusResidual':err,'normalizationBeforeFinalScaling':rawnorm,'elapsedSeconds':time.monotonic()-tic};attempts.append(row)
  if best is None or err<best['maximumReductionFrobeniusResidual']:
   best={**row,'reductions':checks,'coefficients':[[float(z.real),float(z.imag)]for z in v]}
  report={'n':n,'t':t,'dimension':1,'threshold':1e-10,'requestedRestarts':restarts,'maxFunctionEvaluationsPerAttempt':steps,'method':'scipy.optimize.least_squares, trust-region reflective, analytic Jacobian','status':'candidate-found'if best['maximumReductionFrobeniusResidual']<1e-10 else'not-found-in-finite-search','analyticJacobianDirectionalError':jacerr,'seed':seed,'attempts':attempts,'best':best,'elapsedSeconds':time.monotonic()-beg,'scope':'Full-complex finite multistart search. Numerical residuals do not certify exact existence; unsuccessful searches do not prove nonexistence.'}
  (HERE/f'n{n}-t{t}-search.json').write_text(json.dumps(report,indent=2)+'\n')
  print(n,t,r,sol.nfev,err,round(time.monotonic()-beg,2),flush=True)
  if err<2e-12:break
 return report

def main():
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--t',type=int,required=True);p.add_argument('--restarts',type=int,default=24);p.add_argument('--steps',type=int,default=450);p.add_argument('--seed',type=int,default=211103);a=p.parse_args();search(a.n,a.t,a.restarts,a.steps,a.seed+100*a.n+a.t)
if __name__=='__main__':main()
