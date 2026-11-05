"""AC2 search restricted to one Dicke parity.

This is equivalent to searching AC2 children of variable maximal Jz-neutral
parents inside that parity sector; Jx and Jy compressions vanish automatically.
Numerical witnesses still require exact reconstruction.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent

def problem(n,k,parity=0):
    D=n+1;j=n/2;indices=np.arange(parity,D,2);d=len(indices)
    z=np.diag(j-np.arange(D))/j
    plus=np.diag(np.sqrt(np.arange(1,D)*(D-np.arange(1,D))),1)/j
    plus2=plus@plus;c=(j+1)/(3*j)
    full=np.array([np.eye(D),z,z@z-c*np.eye(D),(plus2+plus2.T)/2,(plus2-plus2.T)/(2j)])
    ops=full[:,indices][:,:,indices];ii=np.arange(k);aa,bb=np.triu_indices(k,1);eye=np.eye(k)
    def unpack(v):return v[:d*k].reshape(d,k)+1j*v[d*k:].reshape(d,k)
    def pack(A):return np.concatenate([A[:,ii,ii].real,np.sqrt(2)*A[:,aa,bb].real,np.sqrt(2)*A[:,aa,bb].imag],axis=1).reshape(-1)
    def fun(v):
        B=unpack(v);C=np.einsum('ra,lrs,sb->lab',B.conj(),ops,B);C[0]-=eye;return pack(C)
    def jac(v):
        B=unpack(v);Z=ops@B
        row=np.einsum('ac,lrb->labrc',eye,Z);col=np.einsum('bc,lra->labrc',eye,Z.conj())
        deriv=np.concatenate([(row+col).reshape(5,k,k,-1),(1j*(col-row)).reshape(5,k,k,-1)],axis=-1)
        return np.concatenate([deriv[:,ii,ii].real,np.sqrt(2)*deriv[:,aa,bb].real,np.sqrt(2)*deriv[:,aa,bb].imag],axis=1).reshape(5*k*k,-1)
    return fun,jac,unpack,d,indices

def run(n,k,parity=0,restarts=3,steps=800,seed=5432):
    fun,jac,unpack,d,indices=problem(n,k,parity);rng=np.random.default_rng(seed);best=1e9;start=time.time()
    for attempt in range(restarts):
        B=np.linalg.qr(rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)))[0]
        v=np.r_[B.real.ravel(),B.imag.ravel()]
        if attempt==0:
            direction=rng.normal(size=len(v));h=1e-6
            err=np.max(abs((fun(v+h*direction)-fun(v-h*direction))/(2*h)-jac(v)@direction));assert err<1e-7,err
        sol=least_squares(fun,v,jac=jac,max_nfev=steps,ftol=1e-12,gtol=1e-12,xtol=1e-12)
        residual=float(np.linalg.norm(sol.fun));print(n,k,parity,attempt,sol.nfev,residual,round(time.time()-start,1),flush=True)
        if residual<best:
            best=residual;B=np.zeros((n+1,k),complex);B[indices]=unpack(sol.x)
            data={'n':n,'k':k,'parity':parity,'residual':residual,'attempt':attempt,'seed':seed,'basis':[[[float(z.real),float(z.imag)] for z in row] for row in B],'scope':'Numerical witness only until exact certification. Searches full complex subspaces within a single Dicke parity.'}
            (HERE/f'search-parity-n{n}-k{k}-p{parity}.json').write_text(json.dumps(data,indent=2)+'\n')
        if residual<1e-10:break
    return best

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--parity',type=int,default=0);p.add_argument('--restarts',type=int,default=3);p.add_argument('--steps',type=int,default=800);a=p.parse_args();run(a.n,a.k,a.parity,a.restarts,a.steps)
