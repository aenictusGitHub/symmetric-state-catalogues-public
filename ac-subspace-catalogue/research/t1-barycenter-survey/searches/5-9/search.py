from pathlib import Path
import json,math,itertools,time
import numpy as np
import mpmath as mp
import sympy as sy
from scipy.optimize import minimize,differential_evolution
ROOT=Path(__file__).resolve().parent
SRC=Path(__file__).resolve().parents[3]/'analytical-constructions/t1-bases.json'
mp.mp.dps=90
ROWS={r['n']:r for r in json.loads(SRC.read_text())['subspaces']if 5<=r['n']<=9}
def cstr(z):return [mp.nstr(mp.re(z),80),mp.nstr(mp.im(z),80)]
def xyz(z):
 q=1+abs(z)**2;return [2*z.real/q,-2*z.imag/q,(abs(z)**2-1)/q]
def npbar_roots(z,n,north=0,south=0):
 q=1+abs(z)**2
 return np.array([np.sum(2*z.real/q),np.sum(-2*z.imag/q),np.sum((abs(z)**2-1)/q)+north-south])/n
class Problem:
 def __init__(self,n):
  self.row=ROWS[n];self.n=n;self.k=len(self.row['basis']);self.V=np.array([[complex(*c)for c in b['coefficients']]for b in self.row['basis']]).T
  self.P=self.V*np.array([(-1)**r*math.sqrt(math.comb(n,r))for r in range(n+1)])[:,None]
  self.MV=mp.matrix([[mp.mpc(str(sy.N(sy.re(sy.sympify(self.row['basis'][j]['exactCoefficients'][r])),95)),str(sy.N(sy.im(sy.sympify(self.row['basis'][j]['exactCoefficients'][r])),95)))for j in range(self.k)]for r in range(n+1)])
  self.MP=mp.matrix([[(-1)**r*mp.sqrt(mp.binomial(n,r))*self.MV[r,j]for j in range(self.k)]for r in range(n+1)])
 def bary(self,c):
  p=self.P@c;ix=np.flatnonzero(abs(p)>1e-25*np.linalg.norm(p));lo,hi=ix[0],ix[-1];roots=np.roots(p[lo:hi+1][::-1]);return npbar_roots(roots,self.n,self.n-hi,lo)
 def obj(self,x):
  c=x[:self.k]+1j*x[self.k:];b=self.bary(c);return -np.dot(b,b)
 def ker(self,z,m):
  M=np.array([[sum(self.P[r,j]*math.factorial(r)/math.factorial(r-a)*z**(r-a)for r in range(a,self.n+1))for j in range(self.k)]for a in range(m)])
  M/=np.linalg.norm(M,axis=1)[:,None];_,s,h=np.linalg.svd(M);return h[-1].conj(),s
 def defbar(self,z,m):
  c,s=self.ker(z,m);p=self.P@c;ix=np.flatnonzero(abs(p)>1e-15*np.linalg.norm(p));lo,hi=ix[0],ix[-1]
  q,rem=np.polynomial.polynomial.polydiv(p[:hi+1],np.polynomial.polynomial.polypow([-z,1],m))
  roots=np.r_[np.repeat(z,m),np.roots(q[::-1])];b=npbar_roots(roots,self.n,self.n-hi)
  return np.linalg.norm(b),c,np.linalg.norm(rem),b
 def mpstate(self,z,m,kind='wronskian'):
  M=mp.matrix([[sum(self.MP[r,j]*mp.factorial(r)/mp.factorial(r-a)*z**(r-a)for r in range(a,self.n+1))for j in range(self.k)]for a in range(m)])
  # m=k or k-1; first k-1 constraints uniquely determine the ray.
  A=M[:self.k-1,:];cof=[]
  for j in range(self.k):cof.append((-1)**j*mp.det(mp.matrix([[A[a,b]for b in range(self.k)if b!=j]for a in range(self.k-1)])))
  c=mp.matrix(cof);norm=mp.sqrt(sum(abs(v)**2 for v in c));c=c/norm
  pivot=max(range(self.k),key=lambda j:abs(c[j]));c*=mp.exp(-1j*mp.arg(c[pivot]));d=self.MV*c;p=[(self.MP*c)[r]for r in range(self.n+1)]
  residual=max(abs(v)for v in M*c)
  hi=max(r for r,v in enumerate(p)if abs(v)>mp.mpf('1e-75'));lo=min(r for r,v in enumerate(p)if abs(v)>mp.mpf('1e-75'))
  q=p[:hi+1];remainders=[]
  for _ in range(m):
   b=[q[-1]]
   for v in reversed(q[:-1]):b.append(v+z*b[-1])
   remainders.append(abs(b[-1]));q=list(reversed(b[:-1]))
  lowq=next(i for i,v in enumerate(q)if abs(v)>mp.mpf('1e-65'))
  rsmall=mp.polyroots(list(reversed(q[lowq:])),maxsteps=1000,error=False)
  rr=[z]*m+[mp.mpc(0)]*lowq+list(rsmall);points=[xyz(v)for v in rr]+[[mp.mpf(0),mp.mpf(0),mp.mpf(1)]]*(self.n-hi)
  bar=[sum(v[a]for v in points)/self.n for a in range(3)];distance=mp.sqrt(sum(v*v for v in bar))
  roots=[{'root':cstr(z),'multiplicity':m}]+[{'root':cstr(r),'multiplicity':1}for r in rsmall]
  if lowq:roots.append({'root':['0','0'],'multiplicity':lowq})
  if self.n-hi:roots.append({'root':None,'multiplicity':self.n-hi})
  return dict(n=self.n,dimension=self.k,subspaceId=self.row['id'],distance=mp.nstr(distance,80),barycenter=list(map(lambda x:mp.nstr(x,80),bar)),basisCoordinates=[cstr(v)for v in c],dickeCoefficients=[cstr(v)for v in d],roots=roots,multiplicityPattern=sorted([r['multiplicity']for r in roots],reverse=True),method=kind,derivativeResidual=mp.nstr(residual,12),deflationResidual=mp.nstr(max(remainders),12),optimalityProved=False)
 def wronskian(self):
  k=self.k;n=self.n;out={};support=[[r for r in range(n+1)if self.MP[r,j]!=0]for j in range(k)]
  for ix in itertools.product(*support):
   if len(set(ix))!=k:continue
   coef=mp.mpc(1)
   for j,r in enumerate(ix):coef*=self.MP[r,j]
   for a in range(k):
    for b in range(a+1,k):coef*=ix[b]-ix[a]
   degree=sum(ix)-k*(k-1)//2;out[degree]=out.get(degree,mp.mpc(0))+coef
  scale=max(abs(v)for v in out.values());out={r:c/scale for r,c in out.items()if abs(c)>mp.mpf('1e-80')*scale};lo,hi=min(out),max(out)
  q=[out.get(r,mp.mpc(0))for r in range(lo,hi+1)];roots=np.roots(np.array([complex(z)for z in q])[::-1]);fn=lambda z:mp.polyval(q[::-1],z)
  found=[];bad=[]
  for rr in roots:
   try:
    z=mp.findroot(fn,mp.mpc(rr),tol=mp.mpf('1e-75'),maxsteps=200)
    if any(abs(z-x)<mp.mpf('1e-25')for x in found):continue
    found.append(z)
   except Exception as e:bad.append(str(e))
  return found,{'minDegree':lo,'maxDegree':hi,'finiteNonzeroRoots':len(found),'polynomialDegreeAfterZeroRemoval':hi-lo,'failedRoots':bad,'coefficients':[cstr(v)for v in q]}

def save(name,data): (ROOT/name).write_text(json.dumps(data,indent=2)+'\n')
def run(n):
 P=Problem(n);t=time.time();rng=np.random.default_rng(51720+n);search=[]
 roots,wr=P.wronskian();candidates=[]
 for z in roots:
  try:candidates.append(P.mpstate(z,P.k))
  except Exception as ex:print('wr failure',n,complex(z),str(ex),flush=True)
 save(f'n{n}-wronskian.json',dict(polynomial=wr,candidates=candidates));best=max(candidates,key=lambda r:float(r['distance']),default=None)
 print('WR',n,'best',best['distance']if best else None,flush=True)
 for i in range(240):
  x=rng.normal(size=2*P.k);x/=np.linalg.norm(x);r=minimize(P.obj,x,method='BFGS',options={'maxiter':500,'gtol':1e-10});c=r.x[:P.k]+1j*r.x[P.k:];c/=np.linalg.norm(c);b=P.bary(c)
  search.append(dict(method='BFGS',start=i,distance=float(np.linalg.norm(b)),coordinates=[[z.real,z.imag]for z in c],iterations=r.nit))
 save(f'n{n}-search.json',search);print('BFGS',n,max(x['distance']for x in search),'time',time.time()-t,flush=True)
 for chart in range(P.k):
  def unpack(x):
   c=np.insert(x[:P.k-1]+1j*x[P.k-1:],chart,1+0j);return c/np.linalg.norm(c)
  r=differential_evolution(lambda x:-np.linalg.norm(P.bary(unpack(x))),[(-2,2)]*(2*P.k-2),seed=600+n+chart,maxiter=600,popsize=18,tol=1e-10,polish=False)
  polished=minimize(lambda x:-np.linalg.norm(P.bary(unpack(x))),r.x,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-11,'fatol':1e-13})
  c=unpack(polished.x);search.append(dict(method='DE+NM',chart=chart,distance=-float(polished.fun),coordinates=[[z.real,z.imag]for z in c],evaluations=r.nfev+polished.nfev))
  save(f'n{n}-search.json',search)
 print('DE',n,max(x['distance']for x in search),'time',time.time()-t,flush=True)
 for seed in range(6):
  def obj(x):return -P.defbar(complex(*x),P.k-1)[0]
  rr=differential_evolution(obj,[(-3,3)]*2,seed=800+n*10+seed,maxiter=300,popsize=18,tol=1e-10,polish=True)
  d,c,e,b=P.defbar(complex(*rr.x),P.k-1);search.append(dict(method='DE-repeated',multiplicity=P.k-1,seed=seed,root=rr.x.tolist(),distance=float(d),coordinates=[[z.real,z.imag]for z in c],deflationResidual=float(e),evaluations=rr.nfev))
  save(f'n{n}-search.json',search)
 print('CONSTRAINED',n,max(x['distance']for x in search),'time',time.time()-t,flush=True)
 overall=max(search,key=lambda r:r['distance']);record={'n':n,'subspaceId':P.row['id'],'dimension':P.k,'bestWronskian':best,'bestSearch':overall,'searchStarts':len(search),'timeSeconds':time.time()-t}
 save(f'n{n}-summary.json',record)
if __name__=='__main__':
 import sys
 for n in map(int,sys.argv[1:] or list(range(5,10))):run(n)
