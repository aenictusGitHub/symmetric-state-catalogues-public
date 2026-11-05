from search import *
from scipy.linalg import null_space
for n in range(15,21):
 t=time.time();s=next(x for x in DATA['subspaces']if x['n']==n);k=s['dimension'];V=np.array([[complex(*s['basis'][j]['coefficients'][r])for j in range(k)]for r in range(n+1)])
 P=V*np.array([(-1)**r*math.sqrt(math.comb(n,r))for r in range(n+1)])[:,None];results=[];rng=np.random.default_rng(86000+n)
 def evaluate(c):
  p=P@c;ix=np.flatnonzero(abs(p)>1e-12*np.linalg.norm(p));lo,hi=ix[0],ix[-1];z=np.roots(p[lo:hi+1][::-1]);rad=abs(z)**2;q=1+rad
  b=np.array([sum(2*z.real/q),sum(-2*z.imag/q),sum((rad-1)/q)+n-hi-lo])/n;return b,lo,hi
 seen=set()
 for pole in ['zero','infinity']:
  for m in range(1,n//2+1):
   A=V[:m]if pole=='zero'else V[n+1-m:];U=null_space(A);l=U.shape[1]
   if l==0:continue
   key=(pole,l)
   if key in seen:continue
   seen.add(key)
   def objective(x):
    c=U@(x[:l]+1j*x[l:]);b,_,_=evaluate(c);return -b@b
   count=1 if l==1 else 100
   best=None
   for start in range(count):
    x=rng.normal(size=2*l);res=minimize(objective,x,method='BFGS',options={'gtol':1e-10,'maxiter':600})if l>1 else type('R',(),dict(x=x,nfev=1))()
    c=U@(res.x[:l]+1j*res.x[l:]);c/=np.linalg.norm(c);b,lo,hi=evaluate(c);d=np.linalg.norm(b)
    if best is None or d>best['distance']:best={'pole':pole,'constraintMultiplicity':m,'nullDimension':l,'zeroMultiplicity':int(lo),'infinityMultiplicity':int(n-hi),'distance':float(d),'barycenter':b.tolist(),'basisCoordinates':[[x.real,x.imag]for x in c],'attempts':count}
   results.append(best)
   print(n,pole,m,l,best['distance'],flush=True)
 (ROOT/f'n{n}-poles.json').write_text(json.dumps({'n':n,'results':results,'attempts':sum(r['attempts']for r in results),'best':max(results,key=lambda r:r['distance']),'elapsedSeconds':time.time()-t},indent=2)+'\n')
