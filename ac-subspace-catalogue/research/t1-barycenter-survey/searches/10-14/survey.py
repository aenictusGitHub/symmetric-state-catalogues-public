from pathlib import Path
import json,math,itertools,time,sys
import numpy as np
import mpmath as mp
import sympy as sy
from scipy.optimize import minimize,differential_evolution
ROOT=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve().parents[3]/'analytical-constructions/t1-bases.json'
DATA=json.loads(SOURCE.read_text())
mp.mp.dps=100

def rec(z):return [mp.nstr(mp.re(z),85),mp.nstr(mp.im(z),85)]
def mpval(s):
 x=sy.sympify(s);return mp.mpc(str(sy.N(sy.re(x),105)),str(sy.N(sy.im(x),105)))
def asmp(z):return mp.mpc(float(np.real(z)),float(np.imag(z)))

def get(n):
 s=next(s for s in DATA['subspaces']if s['n']==n);k=s['dimension']
 V=np.array([[complex(*c)for c in v['coefficients']]for v in s['basis']]).T
 M=[[mpval(v['exactCoefficients'][r])for v in s['basis']]for r in range(n+1)]
 P=[[(-1)**r*mp.sqrt(mp.binomial(n,r))*M[r][j]for j in range(k)]for r in range(n+1)]
 return s,V,M,P

def points_bary(roots,n,infinity=0):
 xyz=[[2*z.real/(1+abs(z)**2),-2*z.imag/(1+abs(z)**2),(abs(z)**2-1)/(1+abs(z)**2)]for z in roots]+[[mp.mpf(0),mp.mpf(0),mp.mpf(1)]]*infinity
 b=[sum(row[a]for row in xyz)/n for a in range(3)];return xyz,b,mp.sqrt(sum(x*x for x in b))
def wronskian(P,k,n):
 supports=[[r for r in range(n+1)if P[r][j]!=0]for j in range(k)];W={}
 for rs in itertools.product(*supports):
  v=mp.mpc(1)
  for i in range(k):
   v*=P[rs[i]][i]
   for j in range(i):v*=rs[i]-rs[j]
  ex=sum(rs)-k*(k-1)//2
  if v:W[ex]=W.get(ex,mp.mpc(0))+v
 mx=max(abs(x)for x in W.values());W={e:v for e,v in W.items()if abs(v)>mx*mp.mpf('1e-85')}
 lo=min(W);hi=max(W);w=[W.get(i,mp.mpc(0))/mx for i in range(lo,hi+1)]
 return lo,hi,w

def degenerate_record(n,s,P,Vmp,root,k):
 cols=len(P[0]);D=mp.matrix([[sum(P[r][j]*mp.factorial(r)/mp.factorial(r-l)*root**(r-l)for r in range(l,n+1))for j in range(cols)]for l in range(k-1)])
 dnp=np.array([[complex(D[a,j])for j in range(cols)]for a in range(k-1)]);dnp/=np.linalg.norm(dnp,axis=1)[:,None]
 _,_,vh=np.linalg.svd(dnp);v=vh[-1].conj();pivot=int(np.argmax(abs(v)));keep=[j for j in range(cols)if j!=pivot]
 sol=mp.lu_solve(mp.matrix([[D[a,j]for j in keep]for a in range(k-1)]),mp.matrix([-D[a,pivot]for a in range(k-1)]))
 alpha=[mp.mpc(1)if j==pivot else sol[keep.index(j)]for j in range(cols)]
 norm=mp.sqrt(sum(abs(x)**2 for x in alpha));phase=mp.exp(-1j*mp.arg(alpha[pivot]));alpha=[x*phase/norm for x in alpha]
 coeff=[sum(Vmp[r][j]*alpha[j]for j in range(cols))for r in range(n+1)];poly=[sum(P[r][j]*alpha[j]for j in range(cols))for r in range(n+1)]
 hi=max(r for r in range(n+1)if any(P[r][j]!=0 for j in range(cols)));q=poly[:hi+1];errs=[]
 for _ in range(k):
  b=[q[-1]]
  for a in reversed(q[:-1]):b.append(a+root*b[-1])
  errs.append(abs(b[-1]));q=list(reversed(b[:-1]))
 rr=mp.polyroots(list(reversed(q)),maxsteps=4000,error=False,extraprec=100)
 roots=[root]*k+rr;xyz,b,d=points_bary(roots,n,n-hi)
 A=mp.matrix([[coeff[a+b]*mp.sqrt(mp.binomial(n-1,b)/mp.binomial(n,a+b))for b in range(n)]for a in range(2)])
 rho=A*A.H;ac=max(abs(rho[a,b]-int(a==b)/2)for a in range(2)for b in range(2))
 distinct=[{'root':rec(root),'multiplicity':k}]+[{'root':rec(z),'multiplicity':1}for z in rr]+([{'root':None,'multiplicity':n-hi}]if n-hi else[])
 return dict(n=n,subspaceId=s['id'],dimension=cols,distance=mp.nstr(d,85),barycenter=[mp.nstr(x,85)for x in b],basisCoordinates=[rec(x)for x in alpha],dickeCoefficients=[rec(x)for x in coeff],distinctRoots=distinct,stars=[[mp.nstr(x,85)for x in row]for row in xyz],method='Wronskian root with multiplicity equal to subspace dimension; 100-digit refinement and deflation',deflationResidual=mp.nstr(max(errs),20),oneQubitReductionResidual=mp.nstr(ac,20),optimalityProved=False)

def npbary(alpha,V,P):
 p=P@alpha;n=len(p)-1;ids=np.flatnonzero(abs(p)>np.linalg.norm(p)*1e-28);lo,hi=ids[0],ids[-1]
 z=np.roots(p[lo:hi+1][::-1]);r=abs(z)**2;q=1+r
 return np.array([np.sum(2*z.real/q),np.sum(-2*z.imag/q),np.sum((r-1)/q)+n-hi-lo])/n

def run_wronskian(n):
 s,V,M,P=get(n);k=s['dimension'];lo,hi,w=wronskian(P,k,n)
 zz=np.roots(np.array([complex(x)for x in w[::-1]]));out=[]
 print('N',n,'Wdegree',hi-lo,'zero',lo,flush=True)
 for i,z in enumerate(zz):
  try:
   f=lambda x:mp.polyval(w[::-1],x);root=mp.findroot(f,(asmp(z),asmp(z)*(1+mp.mpf('1e-7'))+mp.mpf('1e-8')),tol=mp.mpf('1e-90'),maxsteps=100)
   if any(abs(root-mp.mpc(*r['root']))<mp.mpf('1e-30')for r in out):continue
   cand=degenerate_record(n,s,P,M,root,k);assert mp.mpf(cand['deflationResidual'])<mp.mpf('1e-65')
   out.append({'root':rec(root),'result':cand});print(n,i,'d',cand['distance'][:24],flush=True)
  except Exception as e:print('FAIL',n,i,type(e).__name__,str(e)[:100],flush=True)
 (ROOT/f'n{n}-wronskian.json').write_text(json.dumps({'n':n,'zeroMultiplicity':lo,'rawDegree':hi,'polynomialCoefficients':[rec(x)for x in w],'candidates':out},indent=2))
 best=max(out,key=lambda r:float(r['result']['distance']))['result'];(ROOT/f'n{n}-result.json').write_text(json.dumps(best,indent=2));return best

def run_search(n):
 s,V,M,Pmp=get(n);k=s['dimension'];P=np.array([[complex(x)for x in row]for row in Pmp]);rng=np.random.default_rng(9031+n)
 def objective(x):
  a=x[:k]+1j*x[k:];b=npbary(a,V,P);return -np.dot(b,b)
 out=[]
 for i in range(220):
  x=rng.normal(size=2*k);x/=np.linalg.norm(x);res=minimize(objective,x,method='BFGS',options={'maxiter':650,'gtol':1e-10});a=res.x[:k]+1j*res.x[k:];a/=np.linalg.norm(a);b=npbary(a,V,P)
  out.append(dict(method='BFGS',start=i,distance=float(np.linalg.norm(b)),basisCoordinates=[[z.real,z.imag]for z in a],evaluations=int(res.nfev)))
  if i%50==0:print('search',n,i,'best',max(r['distance']for r in out),flush=True)
 for chart in range(k):
  def unpack(x):return np.insert(x[:k-1]+1j*x[k-1:],chart,1+0j)
  def obj(x):
   b=npbary(unpack(x),V,P);return -np.dot(b,b)
  res=differential_evolution(obj,[(-2,2)]*(2*k-2),seed=3491+n*10+chart,popsize=15,maxiter=600,tol=1e-9,polish=False)
  r2=minimize(obj,res.x,method='Nelder-Mead',options={'maxiter':6000,'xatol':1e-11,'fatol':1e-13});a=unpack(r2.x);a/=np.linalg.norm(a);b=npbary(a,V,P)
  out.append(dict(method='DE+Nelder-Mead',chart=chart,distance=float(np.linalg.norm(b)),basisCoordinates=[[z.real,z.imag]for z in a],evaluations=int(res.nfev+r2.nfev)))
  print('DE',n,chart,'best',max(r['distance']for r in out),flush=True)
 (ROOT/f'n{n}-search.json').write_text(json.dumps({'attemptCount':len(out),'best':max(out,key=lambda r:r['distance']),'attempts':out},indent=2))
 return max(out,key=lambda r:r['distance'])
if __name__=='__main__':
 mode=sys.argv[1];ns=list(map(int,sys.argv[2:]))
 for n in ns:
  if mode=='w':run_wronskian(n)
  elif mode=='s':run_search(n)
