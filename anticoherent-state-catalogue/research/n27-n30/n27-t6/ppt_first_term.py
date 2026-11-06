"""Validated native first-leg partial-time-reversal dual, arbitrary copy count.

Exact dyadic positive Gram factors, exact rational scale and Arb Wigner6j.
This is future support until a complete certificate invokes and validates it.
"""
import numpy as np
from validated_bootstrap import A,M,Q,sixj_ball

class PPTFirstTerm:
 def __init__(self,root,cert,fam):
  self.n=cert['n'];self.copies=cert['copies'];assert self.copies>=3
  self.parent=fam[self.copies-1]['meta']['source'];self.qblocks={};self.output=[]
  description=cert['pptFirst'];assert description['split']==[1,self.copies-1]
  scale=Q(description['dualScaleRational']);assert scale>0;self.scale=A(scale)
  z=np.load(root/description['factorFile']);Ls=set()
  for o in description['outputs']:
   L=o['L'];assert L not in Ls and L%2==(self.copies*self.n)%2;Ls.add(L)
   valid=[b for b in self.parent if abs(self.n-b['J'])<=L<=self.n+b['J']]
   assert o['parents']==[dict(J=b['J'],p=b['p'])for b in valid]
   offset={};off=0
   for b in valid:offset[b['J']]=(off,off+b['p']);off+=b['p']
   f=z['L'+str(L)];assert f.ndim==2 and f.shape[0]==off and np.all(np.isfinite(f))
   FF=M(f.tolist());q=self.scale*(FF*FF.transpose());self.output.append(dict(L=L,offset=offset,rankColumns=f.shape[1]))
   for si,S in enumerate(valid):
    a,b=offset[S['J']]
    for T in valid[si:]:
     c,d=offset[T['J']];self.qblocks[L,S['J'],T['J']]=M([[q[a+i,c+j]for j in range(T['p'])]for i in range(S['p'])])

 def ambient(self,B,parent):
  n=self.n;J=B['J'];valid=[b for b in parent if abs(n-b['J'])<=J<=n+b['J']];offset={};off=0
  assert [(b['J'],b['p'])for b in parent]==[(b['J'],b['p'])for b in self.parent]
  for b in valid:offset[b['J']]=(off,off+b['p']);off+=b['p']
  H=M(off,off);phase=-1 if(self.copies*n)%2 else 1
  for si,S in enumerate(valid):
   for T in valid[si:]:
    s,t=S['J'],T['J'];block=M(S['p'],T['p'])
    for o in self.output:
     L=o['L'];q=self.qblocks.get((L,s,t))
     if q is None:continue
     coefficient=phase*A((J+1)*(L+1)).sqrt()*sixj_ball(n,s,J,n,t,L)
     if not coefficient.is_zero():block+=q*coefficient
    a,b=offset[s];c,d=offset[t]
    for i in range(S['p']):
     for j in range(T['p']):
      H[a+i,c+j]=block[i,j]
      if s!=t:H[c+j,a+i]=block[i,j]
  return H
