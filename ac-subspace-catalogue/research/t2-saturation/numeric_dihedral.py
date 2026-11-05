"""Exploratory D_k-symmetric AC2 search. Numerical success is not a certificate."""
import json,argparse,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,linprog
HERE=Path(__file__).resolve().parent

def problem(N,k,central_sign=1):
 j=N/2;m=np.arange(N+1)-j;V=j*(j+1)/3
 if (N+1)%k:raise ValueError('Equal-length residue supports needed')
 L=(N+1)//k
 indep=[];cmap={};sign={}
 for a in range((k+1)//2):
  ls=range(L) if a<k//2 else range((L+1)//2)
  for l in ls:
   r=a+k*l
   if a==k//2 and k%2 and central_sign<0 and r==N-r:continue
   ci=len(indep);indep.append((a,r));cmap[r]=ci;sign[r]=1;cmap[N-r]=ci;sign[N-r]=central_sign if a==k//2 and k%2 else 1
 fixed={i for i,(a,r) in enumerate(indep) if r==a}
 re=np.arange(len(indep));im=np.full(len(indep),-1,dtype=int);im[[i for i in range(len(indep)) if i not in fixed]]=np.arange(len(indep),2*len(indep)-len(fixed));dim=2*len(indep)-len(fixed)
 def unpack(x):
  c=x[:len(indep)].astype(complex)
  for ii,vi in enumerate(im):
   if vi>=0:c[ii]+=1j*x[vi]
  return np.array([sign[r]*c[cmap[r]] if r in cmap else 0j for r in range(N+1)])
 groups=[np.arange(a,N+1,k) for a in range(k)];rows=[]
 for inds in groups:
  for w,t in [(np.ones(len(inds)),1),(m[inds]/j,0),((m[inds]**2-V)/j**2,0)]:rows.append(('d',inds,w,t))
 for q in (1,2):
  for a in range(k):
   src=np.arange(a,N+1-q,k);dest=src+q;w=np.ones(len(src))
   for s in range(q):w*=np.sqrt((N-src-s)*(src+s+1))/j
   rows.append(('c',src,dest,w))
   if q==1:rows.append(('c',src,dest,w*(m[src]+.5)/j))
 count=sum(1 if row[0]=='d' else 2 for row in rows)
 def fun(x,jac=False):
  c=unpack(x);f=np.zeros(count);J=np.zeros((count,dim)) if jac else None;ri=0
  for row in rows:
   if row[0]=='d':
    _,rs,w,t=row;f[ri]=np.dot(abs(c[rs])**2,w)-t
    if jac:
     for r,ww in zip(rs,w):
      if r not in cmap:continue
      ci=cmap[r];ss=sign[r];J[ri,re[ci]]+=2*ww*c[r].real*ss
      if im[ci]>=0:J[ri,im[ci]]+=2*ww*c[r].imag*ss
    ri+=1
   else:
    _,src,dest,w=row;val=np.sum(np.conj(c[dest])*c[src]*w);f[ri:ri+2]=[val.real,val.imag]
    if jac:
     for sr,dr,ww in zip(src,dest,w):
      if sr not in cmap or dr not in cmap:continue
      for r,val in [(sr,np.conj(c[dr])*ww),(dr,c[sr]*ww)]:
       ci=cmap[r];ss=sign[r];J[ri,re[ci]]+=ss*val.real;J[ri+1,re[ci]]+=ss*val.imag
       if im[ci]>=0:
        vv=(1j if r==sr else -1j)*val;J[ri,im[ci]]+=ss*vv.real;J[ri+1,im[ci]]+=ss*vv.imag
    ri+=2
  return J if jac else f
 def initial(rng):
  z=np.zeros(len(indep),dtype=complex)
  for a in range((k+1)//2):
   ids=[i for i,(aa,r) in enumerate(indep) if aa==a]
   if a<k//2:
    rs=np.array([indep[i][1] for i in ids]);mult=np.ones(len(ids))
   else:
    rs=np.array([indep[i][1] for i in ids]);mult=np.array([1 if r==N-r else 2 for r in rs])
   if a<k//2:AA=np.array([mult,m[rs]/j*mult,(m[rs]**2-V)/j**2*mult]);bb=[1,0,0]
   else:AA=np.array([mult,(m[rs]**2-V)/j**2*mult]);bb=[1,0]
   sol=linprog(rng.normal(size=len(ids)),A_eq=AA,b_eq=bb,bounds=(0,None),method='highs');pr=sol.x if sol.success else 1/(mult*len(ids));pr=.95*pr+.05/(mult*len(ids))
   vals=np.sqrt(pr)*np.exp(1j*rng.uniform(0,2*np.pi,len(ids)));vals*=np.exp(-1j*np.angle(vals[0]));z[ids]=vals
  x=np.zeros(dim);x[:len(indep)]=z.real
  for ii,vi in enumerate(im):
   if vi>=0:x[vi]=z[ii].imag
  return x
 return fun,initial,unpack

def run(N,k,restarts,seed,max_nfev,central_sign):
 fun,initial,unpack=problem(N,k,central_sign);rng=np.random.default_rng(seed);best=1e9;t=time.time();out=HERE/f'dihedral-n{N}-k{k}-sign{central_sign}.json'
 for rr in range(restarts):
  sol=least_squares(fun,initial(rng),jac=lambda x:fun(x,True),max_nfev=max_nfev,ftol=1e-13,gtol=1e-13,xtol=1e-13);val=np.linalg.norm(sol.fun)
  print(N,k,central_sign,rr,sol.nfev,val,round(time.time()-t,1),flush=True)
  if val<best:
   best=val;c=unpack(sol.x);out.write_text(json.dumps({'N':N,'k':k,'centralSign':central_sign,'residualNorm':float(val),'components':[[float(z.real),float(z.imag)] for z in c],'parameters':sol.x.tolist(),'singularValues':np.linalg.svd(fun(sol.x,True),compute_uv=False).tolist(),'note':'Numerical search, not existence or exclusion proof.'},indent=2)+'\n')
  if best<1e-10:break
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=30);p.add_argument('--seed',type=int,default=569);p.add_argument('--max-nfev',type=int,default=1500);p.add_argument('--central-sign',type=int,default=1);a=p.parse_args();run(a.N,a.k,a.restarts,a.seed,a.max_nfev,a.central_sign)
