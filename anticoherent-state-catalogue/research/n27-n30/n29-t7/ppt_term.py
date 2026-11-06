"""Validated pair-pair recoupling and positive partial-time-reversal dual."""
import numpy as np
from validated_bootstrap import A,M,sixj_ball,frob_upper
class PPTTerm:
 def __init__(self,root,cert,fam):
  assert cert['copies']==4 and cert['n']%2==1;self.n=cert['n'];self.fam=fam;self.pairs=[b['J']for b in fam[2]['meta']['source']];self.outputs=cert['ppt']['outputs'];assert len({o['L']for o in self.outputs})==len(self.outputs);z=np.load(root/cert['ppt']['factorFile']);self.Q={};self.indices={}
  for o in self.outputs:
   L=o['L'];rows=[(P,Q)for P in self.pairs for Q in self.pairs if abs(P-Q)<=L<=P+Q];assert o['labels']==[list(v)for v in rows];v=z['L'+str(L)];assert v.ndim==2 and v.shape[0]==len(rows)and np.all(np.isfinite(v));F=M(v.tolist());self.Q[L]=F*F.transpose();self.indices[L]={row:i for i,row in enumerate(rows)}
 def ambient(self,B,parent):
  n=self.n;J=B['J'];rows=[(P,Q)for P in self.pairs for Q in self.pairs if abs(P-Q)<=J<=P+Q];m=len(rows);K=M(m,m)
  for ai,(P,Q)in enumerate(rows):
   for bi in range(ai,m):
    R,S=rows[bi];v=A(0)
    for L,qq in self.Q.items():
     ix=self.indices[L];i=ix.get((R,Q));j=ix.get((P,S))
     if i is None or j is None:continue
     k=A((J+1)*(L+1)).sqrt()*sixj_ball(P,Q,J,R,S,L)
     if not k.is_zero():v+=k*qq[i,j]
    K[ai,bi]=v;K[bi,ai]=v
  valid=[S for S in parent if abs(n-S['J'])<=J<=n+S['J']];offset={};off=0
  for S in valid:offset[S['J']]=(off,off+S['p']);off+=S['p']
  F=M(m,off);low=self.fam[3]['meta']['lower']
  for ai,(P,Q)in enumerate(rows):
   for S in valid:
    spin=S['J']
    if not abs(n-Q)<=spin<=n+Q:continue
    c=(-1 if(n+(Q+J)//2)%2 else 1)*A((spin+1)*(P+1)).sqrt()*sixj_ball(n,n,P,Q,J,spin)
    if c.is_zero():continue
    e3=self.fam[3]['full'][spin];a,b=offset[spin];qrow=low[str(Q)][0]
    for k in range(S['p']):F[ai,a+k]=c*A(float(e3[qrow,k]))
  return F.transpose()*K*F,frob_upper(K)
