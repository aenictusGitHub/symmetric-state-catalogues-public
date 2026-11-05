"""Numerical rendering of the exact strengthened radical prescription.
Existence for all N is certified separately by integer polynomial inequalities.
"""
from pathlib import Path
import numpy as np,json

def basis(n):
 q,rho=divmod(n,8)
 assert q>=5
 k=q-1;h=2*q+(2 if rho>=6 else 0);offset=2+rho//2 if rho<6 else (rho-2)//2
 j=n/2;x=(j-offset-2*np.arange(k))/h;v=j*(j+1)/(3*h*h)
 w=np.array([v/4+x*x/4-35*x/36+5/6,-v/4-x*x/4+5*x/12+.5,-v/4-x*x/4+13*x/12-.5,v/4+x*x/4-19*x/36+1/6]).T
 if rho==5:w+=(x-1.5)[:,None]/(2*q)*np.array([-1/3,1,-1,1/3])[None,:]
 assert w.min()>0
 R=offset+2*np.arange(k)[:,None]+h*np.arange(4)[None,:]
 phase=np.ones((k,4),complex);margins=[]
 for i in range(k-1):
  r=R[i];L=np.sqrt(w[i]*w[i+1]*(r+1)*(r+2)*(n-r)*(n-r-1))
  a,b,c,d=L;lo=max(abs(a-b),abs(c-d));hi=min(a+b,c+d)
  assert hi>lo
  diagonal=(lo+hi)/2
  def pair(a,b):
   ca=(diagonal*diagonal+a*a-b*b)/(2*diagonal*a);cb=(diagonal*diagonal+b*b-a*a)/(2*diagonal*b)
   return ca+1j*np.sqrt(max(0,1-ca*ca)),cb-1j*np.sqrt(max(0,1-cb*cb))
  z0,z1=pair(a,b);z2,z3=pair(c,d)
  phase[i+1]=phase[i]*np.array([z0,z1,-z2,-z3]);margins.append(float(hi-lo))
 B=np.zeros((n+1,k),complex)
 for i in range(k):B[R[i],i]=np.sqrt(w[i])*phase[i]
 m=j-np.arange(n+1);jp=np.sqrt(np.arange(1,n+1)*(n+1-np.arange(1,n+1)))
 def raise1(X):return np.vstack([jp[:,None]*X[1:],np.zeros((1,X.shape[1]),complex)])
 PB=raise1(B)
 blocks=[B.conj().T@B-np.eye(k),B.conj().T@(m[:,None]*B)/j,B.conj().T@(m[:,None]**2*B)/(j*(j+1)/3)-np.eye(k),B.conj().T@PB/j,B.conj().T@raise1(PB)/(j*j),B.conj().T@(m[:,None]*PB+raise1(m[:,None]*B))/(j*j)]
 residual=max(float(np.max(np.abs(A)))for A in blocks)
 assert residual<1e-11,(n,residual)
 V=np.zeros((n+1,2*k),complex)
 for i in range(k):
  rows=R[i];u=B[rows,i];mm=m[rows];assert min(mm[:2])>0 and max(mm[2:])<0
  pp=np.sqrt(mm[:2])*u[:2];nn=np.sqrt(-mm[2:])*u[2:]
  perp=lambda z:np.array([-z[1].conjugate(),z[0].conjugate()])
  vv=np.r_[perp(pp)/np.sqrt(mm[:2]),perp(nn)/np.sqrt(-mm[2:])]
  vv-=u*np.vdot(u,vv);vv/=np.linalg.norm(vv)
  V[rows,2*i]=u;V[rows,2*i+1]=vv
 parentResidual=max(np.max(np.abs(V.conj().T@V-np.eye(2*k))),np.max(np.abs(V.conj().T@(m[:,None]*V)/j)),np.max(np.abs(V.conj().T@raise1(V)/j)))
 assert parentResidual<1e-11
 return {'n':n,'k':k,'firstOrderParentDimension':2*k,'normalizedParentResidual':float(parentResidual),'normalizedSpinCompressionResidual':residual,'minimumWeight':float(w.min()),'minimumQuadrilateralMargin':min(margins),'supportMaximum':int(R.max())}
if __name__=='__main__':
 rows=[basis(n)for n in list(range(40,161))+[256,512,1000]]
 Path(__file__).with_name('lower-fourband-render-verification.json').write_text(json.dumps({'passed':True,'role':'Supplementary implementation check; all-N proof is exact polynomial certificate','records':rows},indent=2)+'\n')
 print('PASS',len(rows),'sizes; maximum normalized residual',max(r['normalizedSpinCompressionResidual']for r in rows))
