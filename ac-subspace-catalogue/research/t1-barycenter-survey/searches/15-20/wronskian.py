import json,math,itertools,time
from pathlib import Path
import numpy as np
import mpmath as mp
import sympy as sp
mp.mp.dps=100
ROOT=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve().parents[3]/'analytical-constructions/t1-bases.json'
DATA=json.loads(SOURCE.read_text())
def sources(n):return next(s for s in DATA['subspaces']if s['n']==n)
def mpc(c):
 z=sp.sympify(c,rational=True);return mp.mpc(str(sp.N(sp.re(z),100)),str(sp.N(sp.im(z),100)))
def exactV(n):
 s=sources(n);return mp.matrix([[mpc(s['basis'][j]['exactCoefficients'][r])for j in range(s['dimension'])]for r in range(n+1)])
def polyval(p,z):return mp.polyval(list(reversed(p)),z)
def wronskian(n,V):
 k=V.cols;P=mp.matrix([[(-1)**r*mp.sqrt(mp.binomial(n,r))*V[r,j]for j in range(k)]for r in range(n+1)])
 support=[[r for r in range(n+1)if abs(P[r,j])>mp.mpf('1e-90')]for j in range(k)];D=k*(n-k+1);w=[mp.mpc(0)]*(D+1)
 for rs in itertools.product(*support):
  if len(set(rs))<k:continue
  a=mp.mpc(math.prod(rs[j]-rs[i]for i in range(k)for j in range(i+1,k)))
  for j,r in enumerate(rs):a*=P[r,j]
  power=sum(rs)-k*(k-1)//2;w[power]+=a
 scale=max(map(abs,w));w=[x/scale if abs(x)>scale*mp.mpf('1e-80')else mp.mpc(0)for x in w]
 return P,w

def candidate(n,V,P,z,m):
 k=V.cols
 A=mp.matrix([[sum(P[r,j]*mp.factorial(r)/mp.factorial(r-l)*z**(r-l)for r in range(l,n+1))for j in range(k)]for l in range(k-1)])
 normed=np.array(A.tolist(),complex);normed/=np.linalg.norm(normed,axis=1)[:,None]
 free=int(np.argmax(abs(np.linalg.svd(normed)[2][-1])))
 columns=[j for j in range(k)if j!=free];M=mp.matrix([[A[r,j]for j in columns]for r in range(k-1)])
 rest=mp.lu_solve(M,-A[:,free]);c=[mp.mpc(0)]*k;c[free]=mp.mpc(1)
 for j,x in zip(columns,rest):c[j]=x
 norm=mp.sqrt(sum(abs(x)**2 for x in c));phase=mp.exp(-1j*mp.arg(c[free]));c=mp.matrix([x*phase/norm for x in c]);d=V*c;poly=list(P*c)
 support=[r for r,p in enumerate(poly)if abs(p)>mp.mpf('1e-75')];lo,hi=min(support),max(support)
 q=poly[lo:hi+1];remainders=[]
 for _ in range(m):
  b=[q[-1]]
  for x in reversed(q[:-1]):b.append(x+z*b[-1])
  remainders.append(abs(b[-1]));q=list(reversed(b[:-1]))
 # Degree <=13; find ordinary remaining roots from double precision, then refine separately.
 guesses=np.roots(np.array([complex(x)for x in reversed(q)]))
 rr=[]
 for g in guesses:
  r=mp.findroot(lambda x:polyval(q,x),mp.mpc(g),solver='newton',tol=mp.mpf('1e-90'),maxsteps=100)
  rr.append(r)
 assert len(rr)==len(q)-1
 sep=min([abs(a-b)for i,a in enumerate(rr)for b in rr[i+1:]]or[1]);assert sep>mp.mpf('1e-20'),('duplicate remaining root',n,z)
 roots=[z]*m+rr+[mp.mpc(0)]*lo
 points=[[2*x.real/(1+abs(x)**2),-2*x.imag/(1+abs(x)**2),(abs(x)**2-1)/(1+abs(x)**2)]for x in roots]+[[mp.mpf(0),mp.mpf(0),mp.mpf(1)]]*(n-hi)
 assert len(points)==n
 B=[sum(x[i]for x in points)/n for i in range(3)];dist=mp.sqrt(sum(x*x for x in B))
 def st(x):return mp.nstr(x,85)
 def cx(x):return [st(x.real),st(x.imag)]
 # Direct one-qubit reduction.
 split=mp.matrix([[d[a+b]*mp.sqrt(mp.binomial(n-1,b)/mp.binomial(n,a+b))for b in range(n)]for a in range(2)])
 rho=split*split.H;rhoerr=max(abs(rho[a,b]-int(a==b)/2)for a in range(2)for b in range(2))
 record={'n':n,'dimension':k,'subspaceId':sources(n)['id'],'distance':st(dist),'barycenter':list(map(st,B)),'basisCoordinates':list(map(cx,c)),'dickeCoefficients':list(map(cx,d)),'root':cx(z),'roots':[{'root':cx(z),'multiplicity':m}]+[{'root':cx(x),'multiplicity':1}for x in rr]+([{'root':['0','0'],'multiplicity':lo}]if lo else[])+([{'root':None,'multiplicity':n-hi}]if hi<n else[]),'multiplicityPattern':sorted([m]+[1]*len(rr)+([lo]if lo else[])+([n-hi]if hi<n else[]),reverse=True),'method':'Wronskian root, exact source basis evaluated at100 digits; repeated root explicitly deflated.','precisionDigits':100,'globalOptimalityProved':False,'maxDeflationResidual':st(max(remainders)),'oneQubitReductionResidual':st(rhoerr)}
 assert max(remainders)<mp.mpf('1e-65'),record['maxDeflationResidual']
 return record

def run(n):
 t=time.time();V=exactV(n);P,w=wronskian(n,V);s=sources(n);k=s['dimension'];indices=[i for i,x in enumerate(w)if x!=0];low,high=min(indices),max(indices);trim=w[low:high+1]
 guesses=np.roots(np.array([complex(x)for x in reversed(trim)]));candidates=[];fail=[];zs=[]
 for i,g in enumerate(guesses):
  try:
   z=mp.findroot(lambda x:polyval(trim,x),mp.mpc(g),solver='newton',tol=mp.mpf('1e-90'),maxsteps=150)
   if any(abs(z-a)<mp.mpf('1e-50')for a in zs):continue
   zs.append(z);rec=candidate(n,V,P,z,k);candidates.append(rec)
  except Exception as e:fail.append({'guess':[g.real,g.imag],'error':str(e)[:500]})
 candidates.sort(key=lambda r:-float(r['distance']))
 output={'n':n,'subspaceId':s['id'],'dimension':k,'wronskianDegree':high,'wronskianZeroMultiplicity':low,'nonzeroDegree':high-low,'wronskianCoefficients':[[mp.nstr(x.real,85),mp.nstr(x.imag,85)]for x in w],'guesses':len(guesses),'refinedDistinctRoots':len(zs),'successfulCandidates':len(candidates),'failures':fail,'candidates':candidates,'elapsedSeconds':time.time()-t}
 (ROOT/f'n{n}-wronskian.json').write_text(json.dumps(output,indent=2)+'\n')
 if candidates:(ROOT/f'n{n}-candidate.json').write_text(json.dumps(candidates[0],indent=2)+'\n')
 print(n,'best',candidates[0]['distance']if candidates else None,'pattern',candidates[0]['multiplicityPattern']if candidates else None,'candidates',len(candidates),'fail',len(fail),'seconds',time.time()-t,flush=True)
if __name__=='__main__':
 for n in range(15,21):run(n)
