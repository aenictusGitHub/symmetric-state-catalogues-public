from pathlib import Path
import numpy as np,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'t2-general'))
from four_bands import symbolic_weights

def check(n,loss):
 q,rho=divmod(n,8);k=q-loss;h=2*q;offset=loss+1+rho//2
 if k<1:return
 j=n/2;v=j*(j+1)/(3*h*h)
 supports=offset+2*np.arange(k)[:,None]+h*np.arange(4)[None,:]
 x=(j-offset-2*np.arange(k))/h
 W=np.array([[float(w)for w in symbolic_weights(t,v)]for t in x])
 margin=np.min(W)
 if k>1:
  rr=supports[:-1];A=W[:-1]*W[1:]*(rr+1)*(rr+2)*(n-rr)*(n-rr-1)
  if A.min()<0:return n,k,margin,None
  lens=np.sqrt(A);quad=np.min(lens.sum(axis=1)-2*lens.max(axis=1))
 else:quad=0
 return n,k,margin,quad
if __name__=='__main__':
 for loss in [1,2]:
  bad=[]
  for n in range(16,801):
   z=check(n,loss)
   if z and(z[2]<=0 or z[3] is None or z[3]<-1e-8):bad.append(z)
  print('LOSS',loss,'badcount',len(bad),'first',bad[:10],'last',bad[-10:],flush=True)
