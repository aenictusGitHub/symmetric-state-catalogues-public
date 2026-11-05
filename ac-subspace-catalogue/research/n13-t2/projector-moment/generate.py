#!/usr/bin/env python3
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
from itertools import product,permutations
import math,time,json,sys
import numpy as np
from scipy.linalg import null_space,svd
HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('AC_P3_OUTPUT',str(HERE/'generated'))).resolve()
ROOT.mkdir(parents=True,exist_ok=True)
@lru_cache(None)
def cg(a,b,c,ma,mb):
 mc=ma+mb
 if abs(mc)>c or c<abs(a-b)or c>a+b or(a+b+c)%2:return 0.
 j1,j2,J=[Fraction(x,2)for x in(a,b,c)];m1,m2,M=[Fraction(x,2)for x in(ma,mb,mc)]
 vals=[j1+j2-J,j1-j2+J,-j1+j2+J,j1+j2+J+1,j1+m1,j1-m1,j2+m2,j2-m2,J+M,J-M]
 if any(v.denominator!=1 or v<0 for v in vals):return 0.
 A,B,C,D,e,f,g,h,i,l=map(int,vals)
 fac=math.factorial
 rad=Fraction(fac(A)*fac(B)*fac(C),fac(D))*fac(e)*fac(f)*fac(g)*fac(h)*fac(i)*fac(l)*(c+1)
 zz=max(0,int(j2-J-m1),int(j1-J+m2));zu=min(A,f,g)
 value=sum((Fraction((-1)**z,fac(z)*fac(A-z)*fac(f-z)*fac(g-z)*fac(int(J-j2+m1)+z)*fac(int(J-j1-m2)+z))for z in range(zz,zu+1)),Fraction(0))
 # CG phases from 3j cancel: (-1)^(j1-j2+M +j1-j2+M)=+1.
 return float(value)*math.sqrt(float(rad))

def build(n,k,return_basis=False):
 start=time.monotonic();d=n+1
 words={}
 for w in product(range(d),repeat=3):words.setdefault(3*n-2*sum(w),[]).append(w)
 indices={M:{w:i for i,w in enumerate(ws)}for M,ws in words.items()}
 def coupled(J,M):
  aa=[a for a in range(n+1)if abs(2*a-n)<=J<=2*a+n]
  U=np.zeros((len(words[M]),len(aa)))
  for row,(i,j,l)in enumerate(words[M]):
   mi=n-2*i;mj=n-2*j;ml=n-2*l
   for col,a in enumerate(aa):
    U[row,col]=cg(n,n,2*a,mi,mj)*cg(2*a,n,J,mi+mj,ml)
  return aa,U
 full=[];schur=[];vmeta=[]
 for J in sorted(M for M in words if M>=0):
  aa,H=coupled(J,J);mult=len(aa)
  assert np.max(np.abs(H.T@H-np.eye(mult)))<2e-12,(J,'orth')
  A1=np.zeros((mult,mult));A2=np.zeros_like(A1)
  for perm in permutations(range(3)):
   sign=(-1)**sum(perm[a]>perm[b]for a in range(3)for b in range(a+1,3))
   order=[indices[J][tuple(w[a]for a in perm)]for w in words[J]]
   R=H.T@H[order]
   A1+=R/6;A2+=sign*R/6
  def image(A):
   ev,U=np.linalg.eigh((A+A.T)/2);return U[:,ev>.5]
  H1=image(A1);H2=image(A2)
  S12=np.diag([(-1.)**(n-a)for a in aa])
  order23=[indices[J][(w[0],w[2],w[1])]for w in words[J]]
  S23=H.T@H[order23]
  Hplus=image((np.eye(mult)+S12)/2-A1)
  Hminus=(2/np.sqrt(3))*(np.eye(mult)-S12)/2@S23@Hplus
  assert Hminus.shape[1]==0 or np.max(np.abs(Hminus.T@Hminus-np.eye(Hminus.shape[1])))<2e-12
  Uall={M:coupled(J,M)[1]for M in range(-J,J+1,2)}
  ib=len(full);full.append({'J':J,'aa':aa,'U':Uall,'H':H,'S23':S23})
  for name,hs in [('symmetric',[H1]),('antisymmetric',[H2]),('standard',[Hplus,Hminus])]:
   dim=hs[0].shape[1]
   if not dim:continue
   sb=len(schur);info={'J':J,'type':name,'dim':dim,'copies':len(hs),'ib':ib,'H':hs,'variables':[]};schur.append(info)
   for a in range(dim):
    for b in range(a,dim):
     info['variables'].append(len(vmeta));vmeta.append((sb,a,b))
 print('basis',n,'vars',len(vmeta),'schur',[(a['J'],a['type'],a['dim'])for a in schur],'seconds',time.monotonic()-start,flush=True)
 if return_basis:return words,indices,full,schur,vmeta
 # All maps are real. Pair output rows are entries with common total weight.
 pairs={}
 for w in product(range(d),repeat=2):pairs.setdefault(2*n-2*sum(w),[]).append(w)
 pidx={M:{w:i for i,w in enumerate(ws)}for M,ws in pairs.items()}
 rowmeta=[(M,a,b)for M,ws in sorted(pairs.items())for a in range(len(ws))for b in range(a,len(ws))]
 fullmeta=[(ib,a,b)for ib,B in enumerate(full)for a in range(len(B['aa']))for b in range(a,len(B['aa']))]
 nv=len(vmeta);nf=len(fullmeta)
 trace=np.zeros(nv);equations=np.zeros((3*len(rowmeta),nv));theta=np.zeros((nf,nv));direct=np.zeros((nf,nv));iz=np.zeros((3,nf,nv));iztheta=np.zeros((3,nf,nv))
 weight1=np.array([n/2-i for i in range(d)]);weight2=weight1**2-n*(n+2)/12
 pairU={}
 for a in range(n+1):
  M=2*a
  pairU[a]=np.array([cg(n,n,2*a,n-2*i,n-2*j)for i,j in pairs[M]])
 for col,(sb,a,b)in enumerate(vmeta):
  info=schur[sb];B=full[info['ib']];J=B['J'];C=np.zeros((len(B['aa']),len(B['aa'])))
  for hh in info['H']:
   C+=np.outer(hh[:,a],hh[:,b])
   if a!=b:C+=np.outer(hh[:,b],hh[:,a])
  trace[col]=(J+1)*np.trace(C)
  YY={M:U@C@U.T for M,U in B['U'].items()}
  def y(w,z):
   M=3*n-2*sum(w)
   if M not in YY:return 0.
   return YY[M][indices[M][w],indices[M][z]]
  Z={}
  for M,ws in pairs.items():
   A=np.zeros((len(ws),len(ws)))
   for aa0,(j,l)in enumerate(ws):
    for bb0,(q,r)in enumerate(ws):A[aa0,bb0]=sum(y((i,j,l),(i,q,r))for i in range(d))/k
   Z[M]=A
  for rr,(M,aa0,bb0)in enumerate(rowmeta):
   j,l=pairs[M][aa0];q,r=pairs[M][bb0]
   vals=np.array([y((j,i,l),(i,q,r))for i in range(d)])
   equations[rr,col]=sum(vals)-Z[M][aa0,bb0]
   equations[len(rowmeta)+rr,col]=weight1@vals
   equations[2*len(rowmeta)+rr,col]=weight2@vals
  for ib,Bout in enumerate(full):
   M=Bout['J'];ws=words[M];H=Bout['H'];sz=len(ws)
   T=np.zeros((sz,sz));Iz=[np.zeros((sz,sz))for _ in range(3)];It=[np.zeros((sz,sz))for _ in range(3)]
   for rr0,(i,j,l)in enumerate(ws):
    for cc0,(p,q,r)in enumerate(ws):
     T[rr0,cc0]=(-1.)**(i+p)*y((n-p,j,l),(n-i,q,r))
     if i==p:
      m=2*n-2*(j+l);Iz[0][rr0,cc0]=Z[m][pidx[m][(j,l)],pidx[m][(q,r)]]
     if j==q:
      m=2*n-2*(i+l);Iz[1][rr0,cc0]=Z[m][pidx[m][(i,l)],pidx[m][(p,r)]]
     if l==r:
      m=2*n-2*(i+j);Iz[2][rr0,cc0]=Z[m][pidx[m][(i,j)],pidx[m][(p,q)]]
     if i==p:It[0][rr0,cc0]=Iz[0][rr0,cc0]
     if j==q:
      m=2*n-2*(n-p+l);It[1][rr0,cc0]=(-1.)**(i+p)*Z[m][pidx[m][(n-p,l)],pidx[m][(n-i,r)]]
     if l==r:
      m=2*n-2*(n-p+j);It[2][rr0,cc0]=(-1.)**(i+p)*Z[m][pidx[m][(n-p,j)],pidx[m][(n-i,q)]]
   out=H.T@T@H;izout=[H.T@A@H for A in Iz];itout=[H.T@A@H for A in It]
   for fi,(bi,a0,b0)in enumerate(fullmeta):
    if bi!=ib:continue
    theta[fi,col]=out[a0,b0]
    direct[fi,col]=C[a0,b0]if ib==info['ib']else 0.
    for ii in range(3):
     iz[ii,fi,col]=izout[ii][a0,b0];iztheta[ii,fi,col]=itout[ii][a0,b0]
  if col%20==0:print('maps',col,'/',nv,'sec',time.monotonic()-start,flush=True)
 metadata={'n':n,'k':k,'variables':vmeta,'schur':[{key:v for key,v in B.items()if key not in['H']}for B in schur],'full':[{'J':B['J'],'aa':B['aa']}for B in full],'fullmeta':fullmeta,'rowmeta':rowmeta}
 path=ROOT/f'p3-n{n}-k{k}'
 np.savez_compressed(str(path)+'-maps.npz',A=equations,trace=trace,theta=theta,direct=direct,iz=iz,iztheta=iztheta)
 Path(str(path)+'-metadata.json').write_text(json.dumps(metadata))
 print('done',n,nv,'sec',time.monotonic()-start,flush=True)
 return metadata
if __name__=='__main__':build(int(sys.argv[1])if len(sys.argv)>1 else 13,int(sys.argv[2])if len(sys.argv)>2 else 3)
