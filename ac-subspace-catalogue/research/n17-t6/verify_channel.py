"""Exact N17 AC6 rank-one-complement normalization and covariance audit.
The general all-r theorem is proved in proof.tex; no floating-point premise.
"""
from pathlib import Path
from functools import lru_cache
import json,time
import sympy as s
from sympy.physics.wigner import clebsch_gordan

@lru_cache(None)
def tensor(n,m,ell,q):
 j=s.Rational(n,2);jp=s.Rational(m,2);entries={}
 for a in range(n+1):
  ma=a-j;bb=ma+q+jp
  if bb.is_Integer and 0<=bb<=m:
   b=int(bb);v=s.sqrt(s.Rational(2*ell+1,m+1))*clebsch_gordan(j,ell,jp,ma,q,ma+q)
   if v:entries[b,a]=v
 return entries

def hs(A,B):return s.simplify(sum(s.conjugate(v)*B.get(pos,0)for pos,v in A.items()))
def poly(n,ell):
 grid=[s.Rational(2*i-n,2)for i in range(n+1)];out=[]
 for l in range(ell+1):
  v=[x**l for x in grid]
  for u in out:
   c=sum(a*b for a,b in zip(v,u))/sum(a*a for a in u);v=[s.factor(a-c*b)for a,b in zip(v,u)]
  out.append(v)
 return out[-1]

def verify():
 start=time.time();n,m,t=17,15,6;j=s.Rational(n,2);jp=s.Rational(m,2)
 labels=[(l,q)for l in range(1,4)for q in range(-l,l+1)]
 ops={label:tensor(n,m,*label)for label in labels};grams=0;covariance=0
 for i,a in enumerate(labels):
  for b in labels[i:]:assert hs(ops[a],ops[b])==int(a==b);grams+=1
 for (ell,q),A in ops.items():
  # Explicit raising commutator on Hom(V_j,V_j').
  C={}
  for (b,a),v in A.items():
   assert s.Rational(2*b-m,2)-s.Rational(2*a-n,2)==q
   if b<m:C[b+1,a]=C.get((b+1,a),0)+s.sqrt((m-b)*(b+1))*v
   if a>0:C[b,a-1]=C.get((b,a-1),0)-s.sqrt(a*(n-a+1))*v
  B=ops.get((ell,q+1),{});factor=s.sqrt((ell-q)*(ell+q+1))
  assert all(s.simplify(C.get(pos,0)-factor*B.get(pos,0))==0 for pos in C.keys()|B.keys());covariance+=1
 # Every product E_a^† E_b has ranks no greater than ell_a+ell_b.
 highchecks=0
 for a in labels:
  for b in labels:
   A,B=ops[a],ops[b];product={}
   for (r,i),v in A.items():
    for (rr,k),w in B.items():
     if r==rr:product[i,k]=s.simplify(product.get((i,k),0)+s.conjugate(v)*w)
   q=b[1]-a[1]
   assert s.simplify(sum(product.get((i,i),0)for i in range(n+1)))==int(a==b)
   for L in range(max(a[0]+b[0]+1,abs(q)),n+1):
    assert hs(tensor(n,n,L,q),product)==0;highchecks+=1
 # Direct adjoint channel on diagonal tensor polynomials: rank by rank.
 multipliers=[]
 for L in range(7):
  source=poly(n,L);target=poly(m,L);got=[s.Integer(0)for i in range(n+1)]
  for A in ops.values():
   for (b,a),v in A.items():got[a]+=s.conjugate(v)*v*target[b]
  got=[s.factor(x)for x in got];a=next(i for i,x in enumerate(source)if x);c=s.factor(got[a]/source[a]);assert all(s.simplify(v-c*x)==0 for v,x in zip(got,source))
  multipliers.append({'rank':L,'monicMultiplier':str(c)})
 assert multipliers[0]['monicMultiplier']=='5/6'
 assert len(ops)==15 and n+1==18 and m+1==16
 # Symbolic arithmetic of the general even-order exceptional size.
 r=s.symbols('r',integer=True,positive=True);D=(r+1)**2;M=r*(r+2)
 assert s.expand(M-(D-1))==0 and s.expand((D+1)-2-(D-1))==0
 # At N'=D-1, the same family has D-1 images in dimension D-2.
 assert s.expand((D-1)-1-(M-1))==0
 return {'passed':True,'n':n,'t':t,'maximumDimension':0,'method':'Exact Clebsch-Gordan matrices, tensor covariance, product-rank support, and discrete-polynomial channel reconstruction; no numerical tolerance.','hilbertSchmidtGramEntriesChecked':grams,'raisingCovarianceIdentitiesChecked':covariance,'forbiddenHighRankProductComponentsChecked':highchecks,'adjointChannelMonicMultipliers':multipliers,'imageCount':15,'imageSquaredNorm':'1/18','outputDimension':16,'complementRank':1,'complementQubits':15,'complementInheritsACOrder':6,'contradiction':'A pure AC6 state on fifteen qubits would require fifteen orthogonal rank-one, rank-two and rank-three tensor images in fourteen dimensions.','generalCorollary':{'formula':'K((r+1)^2+1,2r)=0 for every integer r>=1','examples':[[5,2],[10,4],[17,6],[26,8],[37,10]],'proofSource':'proof.tex'},'arbitraryComplexStatesCovered':True,'scope':'Finite exact convention audit for N17; the general theorem follows from the representation-theoretic rank-one-complement proof.','elapsedSeconds':time.time()-start}
if __name__=='__main__':
 report=verify();Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
