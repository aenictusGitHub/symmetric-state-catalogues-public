from pathlib import Path
import numpy as np, scipy.linalg as la, json, argparse,time

def setup(n,a):
 j=n/2;d=n+1;Z=np.diag(j-np.arange(d));P=np.diag(np.sqrt(np.arange(1,d)*(n-np.arange(1,d)+1)),1);X=(P+P.T)/2;Y=(P-P.T)/2j
 raw=[X,Y,Z,X@X-Y@Y,3*Z@Z-j*(j+1)*np.eye(d),X@Y+Y@X,X@Z+Z@X,Y@Z+Z@Y]
 T=[]
 for A in raw:
  for B in T:A=A-np.trace(B.conj().T@A).real*B
  T.append(A/np.linalg.norm(A))
 T=np.array(T)
 def project(W):
  w=W.reshape(d,a,d,a);coeff=np.einsum('lxy,xiyj->lij',T.conj(),w,optimize=True)
  return np.einsum('lxy,lij->xiyj',T,coeff,optimize=True).reshape(d*a,d*a)
 s=(a-1)/2;sz=np.diag(s-np.arange(a));sp=np.diag(np.sqrt(np.arange(1,a)*(a-np.arange(1,a))),1);sx=(sp+sp.T)/2;sy=(sp-sp.T)/2j
 Q=np.kron(X,sx)+np.kron(Y,sy)+np.kron(Z,sz)
 # Choose the best scalar quadratic as a starting point.
 const=j*(j+1)*s*(s+1)/3
 ls=np.arange(abs(j-s),j+s+1);xs=(ls*(ls+1)-j*(j+1)-s*(s+1))/2
 cuts=sorted(set((const-x*x)/x for x in xs if abs(x)>1e-10));trials=[cuts[0]-1]+[(u+v)/2 for u,v in zip(cuts,cuts[1:])]+[cuts[-1]+1]
 starts=[]
 for v in trials:
  eig=xs*xs+v*xs-const;pos=int(sum(2*l+1 for l,e in zip(ls,eig)if e>0));neg=d*a-pos
  starts.append((min(pos,neg),v,1 if pos<=neg else -1))
 _,v,sign=min(starts)
 W=sign*(Q@Q+v*Q-const*np.eye(d*a))
 return T,project,W

def main():
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=35);p.add_argument('--a',type=int,default=5);p.add_argument('--steps',type=int,default=200);p.add_argument('--trials',type=int,default=6);args=p.parse_args();n,a=args.n,args.a;D=(n+1)*a;target=(min((n-2)//4,(n+5)//5))*a-1;T,project,W0=setup(n,a);rng=np.random.default_rng(861);best=999;start=time.time()
 for trial in range(args.trials):
  W=W0/np.linalg.norm(W0)*np.sqrt(D)
  if trial:
   h=rng.normal(size=(D,D))+1j*rng.normal(size=(D,D));h=(h+h.conj().T)/2;h=project(h);W+=h/np.linalg.norm(h)*(.2 if trial<3 else 1)*np.sqrt(D)
  for step in range(args.steps):
   e,U=la.eigh(W,check_finite=False,driver='evr');boundary=e[D-target-1]
   if boundary<best:best=boundary
   if boundary< -1e-7:
    path=Path(__file__).with_name(f'inertia-n{n}-a{a}.npz');np.savez(path,W=W,T=T);print('SUCCESS',n,a,target,boundary,flush=True);return
   ep=e.copy();ep[:D-target]=np.minimum(ep[:D-target],-.001);ep[D-target:]=np.maximum(ep[D-target:],.001)
   V=(U*ep)@U.conj().T;W=project(V);W=W/np.linalg.norm(W)*np.sqrt(D)
  print(n,a,'trial',trial,'bestboundary',best,'last',boundary,'seconds',time.time()-start,flush=True)
 print('NO WITNESS',flush=True)
if __name__=='__main__':main()
