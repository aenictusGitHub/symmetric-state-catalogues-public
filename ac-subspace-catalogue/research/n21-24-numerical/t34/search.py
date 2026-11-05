#!/usr/bin/env python3
"""Reproducible complex AC3/AC4 subspace searches; failed searches are not exclusions."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
from pathlib import Path
import json,time,math,argparse
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent;CAT=HERE.parents[2]
def tensors(n,t):
 d=n+1;j=n/2;plus=np.diag(np.sqrt(np.arange(1,d)*(d-np.arange(1,d))),1).astype(complex);minus=plus.conj().T
 ops=[np.eye(d,dtype=complex)]
 for ell in range(1,t+1):
  T=np.linalg.matrix_power(plus,ell);T/=np.linalg.norm(T)
  for q in range(ell,-1,-1):
   if q:
    ops.extend([(T+T.conj().T)*np.sqrt(d/2),(T-T.conj().T)*np.sqrt(d/2)/1j])
    T=(minus@T-T@minus)/np.sqrt((ell+q)*(ell-q+1))
   else:
    assert np.linalg.norm(T-T.conj().T)<1e-12
    ops.append(T*np.sqrt(d))
 ops=np.array(ops);G=np.einsum('aij,bij->ab',ops.conj(),ops)
 assert np.max(abs(G-d*np.eye(len(ops))))<1e-10
 return ops

def problem(n,t,k,mask=None):
 d=n+1;ops=tensors(n,t);num=len(ops);eye=np.eye(k);ii=np.arange(k);aa,bb=np.triu_indices(k,1)
 if mask is None:mask=np.ones((d,k),bool)
 ix=np.flatnonzero(mask.ravel());size=len(ix)
 def unpack(v):
  B=np.zeros(d*k,complex);B[ix]=v[:size]+1j*v[size:];return B.reshape(d,k)
 def coords(B):return np.r_[B.ravel()[ix].real,B.ravel()[ix].imag]
 def pack(C):return np.concatenate([C[:,ii,ii].real,np.sqrt(2)*C[:,aa,bb].real,np.sqrt(2)*C[:,aa,bb].imag],axis=1).reshape(-1)
 def fun(v):
  B=unpack(v);C=np.einsum('ra,lrb->lab',B.conj(),ops@B);C[0]-=eye;return pack(C)
 def jac(v):
  B=unpack(v);Z=ops@B
  row=np.einsum('ac,lrb->labrc',eye,Z);col=np.einsum('bc,lra->labrc',eye,Z.conj())
  deriv=np.concatenate([(row+col).reshape(num,k,k,-1)[...,ix],(1j*(col-row)).reshape(num,k,k,-1)[...,ix]],axis=-1)
  return np.concatenate([deriv[:,ii,ii].real,np.sqrt(2)*deriv[:,aa,bb].real,np.sqrt(2)*deriv[:,aa,bb].imag],axis=1).reshape(num*k*k,-1)
 return fun,jac,unpack,coords

def independent_check(B,n,t):
 Q=np.linalg.qr(B)[0];d=t+1;mats=[np.array([[Q[a+h,i]*math.sqrt(math.comb(t,a)*math.comb(n-t,h)/math.comb(n,a+h))for h in range(n-t+1)]for a in range(d)])for i in range(Q.shape[1])]
 E=[float(np.linalg.norm(A@D.conj().T-(np.eye(d)/d if i==j else 0)))for i,A in enumerate(mats)for j,D in enumerate(mats)]
 return Q,{'orthonormalityFrobenius':float(np.linalg.norm(Q.conj().T@Q-np.eye(Q.shape[1]))),'maxCrossReductionFrobenius':max(E),'allOrderedCrossReductions':E,'orderedComplexEntriesChecked':d*d*Q.shape[1]**2}

def records():
 data=json.loads((CAT/'data/catalogue.json').read_text());return data['subspaces']

def continued_seeds(n,t,k):
 out=[]
 # Reinterpret nearby Dicke arrays after zero padding or one-qubit contractions.
 sources=records()
 for file in sorted(HERE.glob('n*-t*-k*.json')):
  data=json.loads(file.read_text())
  if data.get('found') and (data['n'],data['t'],data['dimension'])!=(n,t,k):
   sources.append({'id':'local-'+file.stem,'n':data['n'],'t':data['t'],'dimension':data['dimension'],'basis':[{'coefficients':row}for row in data['bestBasis']]})
 for r in sources:
  if r['dimension']<k or r['t']<min(3,t) or abs(r['n']-n)>5:continue
  if r['dimension']!=k:continue
  B=np.array([[complex(*z)for z in row['coefficients']]for row in r['basis']]).T
  if r['n']<=n:
   for shift in sorted(set([0,(n-r['n'])//2,n-r['n']])):
    Z=np.zeros((n+1,k),complex);Z[shift:shift+len(B)]=B;out.append((r['id']+f'-pad{shift}',Z))
  else:
   for which in [0,1]:
    Z=B.copy()
    for m in range(r['n'],n,-1):
     Z=(np.sqrt((m-np.arange(m))/m)[:,None]*Z[:-1])if which==0 else(np.sqrt((np.arange(m)+1)/m)[:,None]*Z[1:])
    if np.linalg.matrix_rank(Z)==k:out.append((r['id']+f'-contract{which}',np.linalg.qr(Z)[0]))
 return out

def run(n,t,k,restarts=12,steps=1500,seconds=300,seed=0,out=None):
 folder=Path(out or HERE);folder.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed or 100000*n+1000*t+10*k)
 fun,jac,unpack,coords=problem(n,t,k);B=np.linalg.qr(rng.normal(size=(n+1,k))+1j*rng.normal(size=(n+1,k)))[0];x=coords(B);v=rng.normal(size=len(x));eps=1e-6
 err=float(np.max(abs((fun(x+eps*v)-fun(x-eps*v))/(2*eps)-jac(x)@v)));assert err<1e-7
 seeds=[('random-0',B)]+continued_seeds(n,t,k)
 for i in range(1,restarts):seeds.append((f'random-{i}',np.linalg.qr(rng.normal(size=(n+1,k))+1j*rng.normal(size=(n+1,k)))[0]))
 oldpath=folder/f'n{n}-t{t}-k{k}.json'
 if oldpath.exists():
  old=json.loads(oldpath.read_text());seeds.insert(0,('previous-best',np.array([[complex(*z)for z in row]for row in old['bestBasis']]).T))
 best=np.inf;start=time.monotonic();attempts=[];winner=None
 class Deadline(Exception):pass
 for name,B in seeds:
  if time.monotonic()-start>seconds:break
  local=[float(np.linalg.norm(fun(coords(B)))),coords(B).copy()];calls=[0]
  def f(x):
   if time.monotonic()-start>seconds:raise Deadline()
   value=fun(x);calls[0]+=1;norm=float(np.linalg.norm(value))
   if norm<local[0]:local[:]=[norm,x.copy()]
   return value
  begin=time.monotonic()
  try:
   sol=least_squares(f,coords(B),jac=jac,max_nfev=steps,ftol=1e-12,xtol=1e-12,gtol=1e-12,method='trf');status=int(sol.status)
  except Deadline:status='time-limit'
  Q,check=independent_check(unpack(local[1]),n,t)
  attempt={'seedName':name,'residualL2':local[0],'nfev':calls[0],'status':status,'seconds':time.monotonic()-begin,**check};attempts.append(attempt)
  print(n,t,k,name,calls[0],local[0],check['maxCrossReductionFrobenius'],round(time.monotonic()-start,2),flush=True)
  if check['maxCrossReductionFrobenius']<best:
   best=check['maxCrossReductionFrobenius'];winner=Q;bestname=name
  report={'n':n,'t':t,'dimension':k,'method':'Unrestricted complex orthonormal-subspace least squares with an analytic Jacobian; independent normalized-Dicke reduced-matrix check after QR.','seed':seed or 100000*n+1000*t+10*k,'jacobianDirectionalError':err,'attempts':attempts,'elapsedSeconds':time.monotonic()-start,'bestSeed':bestname,'maxCrossReductionFrobenius':best,'found':best<1e-10,'bestBasis':[[[float(z.real),float(z.imag)]for z in winner[:,i]]for i in range(k)],'scope':'Numerical candidate only. A failed finite search is not a proof of nonexistence.'}
  oldpath.write_text(json.dumps(report,indent=2)+'\n')
  if best<1e-12:break
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--t',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=12);p.add_argument('--steps',type=int,default=1500);p.add_argument('--seconds',type=float,default=300);p.add_argument('--seed',type=int,default=0);p.add_argument('--out');a=p.parse_args();run(**vars(a))
