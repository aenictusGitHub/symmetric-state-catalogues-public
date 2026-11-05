import itertools,json,math,sys,time
from fractions import Fraction
from optimize import *

def sparse(n,k,mode,maxsupport=3):
 obj=Sector(n,k);sgn=-1 if mode=='max' else 1;rng=np.random.default_rng(8300+n*100+k)
 target=json.loads((ROOT/f'n{n}-s{k}-{mode}.json').read_text())['value']
 matches=[];best=None;start=time.time()
 for size in range(1,maxsupport+1):
  supports=list(itertools.combinations(range(n+1),size))
  for support in supports:
   if support[0]>n-support[-1]:continue # Reflection equivalent.
   idx=np.array(list(support)+[i+n+1 for i in support]);base=np.zeros(2*(n+1))
   def f(z):
    x=base.copy();x[idx]=z
    val,grad=obj.value_grad(x);return sgn*val,sgn*grad[idx]
   for j in range(2 if size>1 else 1):
    x=rng.normal(size=size*2)
    res=minimize(f,x,jac=True,method='L-BFGS-B',options={'maxiter':600,'ftol':2e-15,'gtol':1e-8,'maxls':40})
    val=sgn*res.fun
    z=res.x/np.linalg.norm(res.x);c=np.zeros(n+1,complex);c[list(support)]=z[:size]+1j*z[size:];c*=np.exp(-1j*np.angle(c[np.argmax(abs(c))]))
    result={'n':n,'k':k,'mode':mode,'value':float(val),'support':list(support),'coefficients':[[float(a.real),float(a.imag)] for a in c],'weights':[float(abs(a)**2) for a in c],'weightRationals':[str(Fraction(float(abs(a)**2)).limit_denominator(10000)) for a in c],'phases':[float(np.angle(a)) for a in c]}
    if best is None or sgn*val<sgn*best['value']:best=result
    if abs(val-target)<1e-6:
     matches.append(result)
     break
  if matches:break
 (ROOT/f'n{n}-s{k}-{mode}-sparse.json').write_text(json.dumps({'best':best,'matches':matches,'elapsed':time.time()-start},indent=2)+'\n')
 print(f'sparse n{n} S{k} {mode} best={best["value"]:.12g} target={target:.12g} support={best["support"]} matches={len(matches)} time={time.time()-start:.1f}',flush=True)
 return matches

if __name__=='__main__':
 if len(sys.argv)>1:sparse(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],int(sys.argv[4]) if len(sys.argv)>4 else 3)
 else:
  tasks=[(9,6,'min'),(10,8,'max'),(10,9,'max')]+[(11,k,'max') for k in range(6,11)]+[(12,k,'max') for k in range(6,12)]
  for t in tasks:sparse(*t)
