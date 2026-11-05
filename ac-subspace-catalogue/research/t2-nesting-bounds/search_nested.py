"""Exploratory AC2 search in phase-deformed reflected AC1 parents.
Numerical residuals do not establish nonexistence.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent

def problem(n,k):
    assert n%2==0 or n%6==5
    j=n/2;D=n+1
    central=n%2==0 and int(j)%3!=2
    xs=np.arange(2 if central else 1,int(j),3) if n%2==0 else np.arange(1.5,j,3)
    d=2*len(xs)+central
    V0=np.zeros((D,d),complex);masks=[];col=0
    if central:V0[int(j),0]=1;col+=1
    for x in xs:
        a=np.sqrt((x+1)/(2*x+1));b=np.sqrt(x/(2*x+1))
        V0[int(j-x),col]=a;V0[int(j+x+1),col]=b
        V0[int(j-x-1),col+1]=b;V0[int(j+x),col+1]=-a
        mask=np.zeros((D,d));mask[int(j+x+1),col]=1;mask[int(j+x),col+1]=1
        masks.append(mask);col+=2
    masks=np.array(masks)
    z=np.diag(j-np.arange(D))/j
    plus=np.diag(np.sqrt(np.arange(1,D)*(D-np.arange(1,D))),1)/j
    x=(plus+plus.T)/2;y=(plus-plus.T)/(2j);c=(j+1)/(3*j)
    O=np.array([np.eye(D),x@x-y@y,x@y+y@x,x@z+z@x,y@z+z@y,z@z-c*np.eye(D)])
    ii=np.arange(k);aa,bb=np.triu_indices(k,1);eye=np.eye(k)
    def unpack(v):
        X=v[:d*k].reshape(d,k)+1j*v[d*k:2*d*k].reshape(d,k)
        phase=np.einsum('a,ars->rs',v[2*d*k:],masks)
        V=V0*np.exp(1j*phase)
        return X,V
    def pack(A):return np.concatenate([A[:,ii,ii].real,np.sqrt(2)*A[:,aa,bb].real,np.sqrt(2)*A[:,aa,bb].imag],axis=1).reshape(-1)
    def fun(v):
        X,V=unpack(v);B=V@X;C=np.einsum('ra,lrs,sb->lab',B.conj(),O,B);C[0]-=eye
        return pack(C)
    def jac(v):
        X,V=unpack(v);B=V@X;ops=V.conj().T@O@V;Z=ops@X
        row=np.einsum('ac,lrb->labrc',eye,Z);col=np.einsum('bc,lra->labrc',eye,Z.conj())
        deriv=np.concatenate([(row+col).reshape(6,k,k,-1),(1j*(col-row)).reshape(6,k,k,-1)],axis=-1)
        dB=1j*np.einsum('prs,rs,sk->prk',masks,V,X)
        phases=np.einsum('pra,lrs,sb->labp',dB.conj(),O,B)+np.einsum('ra,lrs,psb->labp',B.conj(),O,dB)
        deriv=np.concatenate([deriv,phases],axis=-1)
        return np.concatenate([deriv[:,ii,ii].real,np.sqrt(2)*deriv[:,aa,bb].real,np.sqrt(2)*deriv[:,aa,bb].imag],axis=1).reshape(6*k*k,-1)
    return fun,jac,unpack,d,len(xs)

def run(n,k,restarts=4,steps=700,seed=1234):
    fun,jac,unpack,d,q=problem(n,k);rng=np.random.default_rng(seed);best=1e9;start=time.time()
    for attempt in range(restarts):
        X=np.linalg.qr(rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)))[0]
        v=np.r_[X.real.ravel(),X.imag.ravel(),rng.uniform(-np.pi,np.pi,q)]
        if attempt==0:
            direction=rng.normal(size=len(v));h=1e-6
            err=np.max(abs((fun(v+h*direction)-fun(v-h*direction))/(2*h)-jac(v)@direction));assert err<1e-7,err
        sol=least_squares(fun,v,jac=jac,max_nfev=steps,ftol=1e-11,gtol=1e-11,xtol=1e-11)
        residual=float(np.linalg.norm(sol.fun));print(n,k,attempt,sol.nfev,residual,round(time.time()-start,1),flush=True)
        if residual<best:
            best=residual;X,V=unpack(sol.x);B=V@X
            def complex_list(A):return [[[float(z.real),float(z.imag)]for z in row]for row in A]
            data={'n':n,'k':k,'parentDimension':d,'residual':residual,'attempt':attempt,'seed':seed,'phases':sol.x[2*d*k:].tolist(),'childCoordinates':complex_list(X),'basis':complex_list(B),'scope':'Numerical search in phase-deformed canonical reflection parents; not an exact existence or nonexistence certificate.'}
            (HERE/f'search-n{n}-k{k}.json').write_text(json.dumps(data,indent=2)+'\n')
        if residual<1e-10:break
    return best

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--restarts',type=int,default=4);p.add_argument('--steps',type=int,default=700);a=p.parse_args();run(a.n,a.k,a.restarts,a.steps)
