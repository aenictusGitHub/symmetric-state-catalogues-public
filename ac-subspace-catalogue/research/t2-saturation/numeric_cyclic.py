"""Exploratory cyclic-residue AC2 search; numerical output is not a certificate."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,linprog
HERE=Path(__file__).resolve().parent

def problem(N,k):
    j=N/2;m=np.arange(N+1)-j;c=j*(j+1)/3
    # Independent phase of each codeword makes its first coefficient real.
    xr=np.arange(N+1); yi=np.full(N+1,-1,dtype=int);yi[k:]=np.arange(N+1,2*(N+1)-k)
    dim=2*(N+1)-k
    def unpack(x):
        z=x[:N+1].astype(complex);z[k:]+=1j*x[N+1:];return z
    groups=[np.arange(a,N+1,k) for a in range(k)]
    rows=[]
    for a,inds in enumerate(groups):
        for ws,target in ((np.ones(len(inds)),1),(m[inds]/j,0),((m[inds]**2-c)/j**2,0)):
            rows.append(('diag',inds,ws,target))
    for q in (1,2):
        for a in range(k):
            src=np.arange(a,N+1-q,k);dest=src+q
            w=np.ones(len(src))
            for s in range(q):w*=np.sqrt((N-src-s)*(src+s+1))/j
            rows.append(('cross',src,dest,w))
            if q==1:rows.append(('cross',src,dest,w*(m[src]+.5)/j))
    count=sum(1 if row[0]=='diag' else 2 for row in rows)
    def both(x,jac=False):
        z=unpack(x);f=np.zeros(count);A=np.zeros((count,dim)) if jac else None;ri=0
        for row in rows:
            if row[0]=='diag':
                _,inds,w,target=row;f[ri]=np.sum(abs(z[inds])**2*w)-target
                if jac:
                    A[ri,xr[inds]]=2*w*z[inds].real
                    for idx,ww in zip(inds,w):
                        if yi[idx]>=0:A[ri,yi[idx]]=2*ww*z[idx].imag
                ri+=1
            else:
                _,src,dest,w=row;v=np.sum(np.conj(z[dest])*z[src]*w);f[ri:ri+2]=[v.real,v.imag]
                if jac:
                    for s,d,ww in zip(src,dest,w):
                        for idx,val in ((xr[s],np.conj(z[d])*ww),(xr[d],z[s]*ww),(yi[s],1j*np.conj(z[d])*ww),(yi[d],-1j*z[s]*ww)):
                            if idx>=0:A[ri,idx]+=val.real;A[ri+1,idx]+=val.imag
                ri+=2
        return A if jac else f
    def initial(rng):
        z=np.zeros(N+1,dtype=complex)
        for inds in groups:
            lin=np.stack([np.ones(len(inds)),m[inds]/j,(m[inds]**2-c)/j**2])
            # Minimize randomized linear objective, mix with a uniform positive distribution.
            p=linprog(rng.normal(size=len(inds)),A_eq=lin,b_eq=[1,0,0],bounds=(0,None),method='highs')
            prob=p.x if p.success else np.ones(len(inds))/len(inds)
            prob=.95*prob+.05/len(inds)
            zs=np.sqrt(prob)*np.exp(1j*rng.uniform(0,2*np.pi,len(inds)));zs*=np.exp(-1j*np.angle(zs[0]));z[inds]=zs
        return np.r_[z.real,z[k:].imag]
    return both,initial,unpack

def run(N,k,restarts,seed,max_nfev):
    fun,initial,unpack=problem(N,k);rng=np.random.default_rng(seed);best=1e100;started=time.time()
    out=HERE/f'cyclic-n{N}-k{k}.json'
    for r in range(restarts):
        x=initial(rng)
        sol=least_squares(fun,x,jac=lambda x:fun(x,True),max_nfev=max_nfev,ftol=1e-13,gtol=1e-13,xtol=1e-13)
        val=np.linalg.norm(sol.fun)
        print(N,k,r,sol.nfev,val,'secs',round(time.time()-started,1),flush=True)
        if val<best:
            best=val;z=unpack(sol.x)
            result={'N':N,'k':k,'residualNorm':val,'maxResidual':float(max(abs(sol.fun))),'components':[[float(a.real),float(a.imag)] for a in z],'squaredAmplitudes':[float(abs(a)**2) for a in z],'singularValues':np.linalg.svd(fun(sol.x,True),compute_uv=False).tolist(),'attempt':r,'seed':seed,'note':'Uncertified numerical search; not a proof of existence or exclusion.'}
            out.write_text(json.dumps(result,indent=2)+'\n')
        if val<1e-10:break
    return best

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=8);p.add_argument('--seed',type=int,default=341);p.add_argument('--max-nfev',type=int,default=1500);a=p.parse_args();run(a.N,a.k,a.restarts,a.seed,a.max_nfev)
