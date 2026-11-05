import json,math,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,differential_evolution
ROOT=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve().parents[3]/'analytical-constructions/t1-bases.json'
DATA=json.loads(SOURCE.read_text())
def run(n):
 s=next(x for x in DATA['subspaces']if x['n']==n);k=s['dimension'];V=np.array([[complex(*s['basis'][j]['coefficients'][r])for j in range(k)]for r in range(n+1)])
 P=V*np.array([(-1)**r*math.sqrt(math.comb(n,r))for r in range(n+1)])[:,None]
 def bary(c):
  p=P@c;support=np.flatnonzero(abs(p)>1e-28*np.linalg.norm(p));lo,hi=support[0],support[-1];z=np.roots(p[lo:hi+1][::-1]);rad=abs(z)**2;q=1+rad
  return np.array([np.sum(2*z.real/q),np.sum(-2*z.imag/q),np.sum((rad-1)/q)+n-hi-lo])/n
 def objective(x):
  b=bary(x[:k]+1j*x[k:]);return -np.dot(b,b)
 def row(res,start,method):
  c=res.x[:k]+1j*res.x[k:];c/=np.linalg.norm(c);c*=np.exp(-1j*np.angle(c[np.argmax(abs(c))]));b=bary(c)
  return {'start':start,'method':method,'distance':float(np.linalg.norm(b)),'basisCoordinates':[[x.real,x.imag]for x in c],'barycenter':b.tolist(),'iterations':int(res.nit),'evaluations':int(res.nfev),'success':bool(res.success)}
 rng=np.random.default_rng(481000+n);out=[];t=time.time();best=0
 for a in range(220):
  x=rng.normal(size=2*k);x/=np.linalg.norm(x)
  opt=minimize(objective,x,method='BFGS',options={'gtol':1e-9,'maxiter':500});rec=row(opt,a,'BFGS');out.append(rec)
  if rec['distance']>best+1e-7:best=rec['distance'];print(n,'BFGS',a,best,flush=True)
 (ROOT/f'n{n}-bfgs.json').write_text(json.dumps(out,indent=2)+'\n')
 bestrow=max(out,key=lambda r:r['distance']);c=np.array([complex(*x)for x in bestrow['basisCoordinates']]);chart=int(np.argmax(abs(c)))
 def unpack(x):return np.insert(x[:k-1]+1j*x[k-1:],chart,1.+0j)
 def de_obj(x):
  b=bary(unpack(x));return -np.dot(b,b)
 opt=differential_evolution(de_obj,[(-2,2)]*(2*k-2),seed=52200+n,popsize=14,maxiter=700,tol=1e-9,polish=False)
 finish=minimize(de_obj,opt.x,method='Nelder-Mead',options={'maxiter':12000,'xatol':1e-11,'fatol':1e-13})
 cc=unpack(finish.x);xx=np.r_[cc.real,cc.imag];finish.x=xx;rec=row(finish,0,'Differential evolution then Nelder-Mead');rec.update(deEvaluations=int(opt.nfev),deIterations=int(opt.nit),chart=chart);out.append(rec)
 combined={'n':n,'subspaceId':s['id'],'dimension':k,'attempts':221,'results':out,'best':max(out,key=lambda r:r['distance']),'elapsedSeconds':time.time()-t}
 (ROOT/f'n{n}-search.json').write_text(json.dumps(combined,indent=2)+'\n')
 print(n,'FINISHED',combined['best']['distance'],'seconds',time.time()-t,flush=True)
if __name__=='__main__':
 for n in range(15,21):run(n)
