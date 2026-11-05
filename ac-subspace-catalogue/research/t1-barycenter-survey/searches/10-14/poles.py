from survey import *

def run(n):
 s,V,M,Pmp=get(n);P=np.array([[complex(x)for x in row]for row in Pmp]);k=s['dimension'];groups=[]
 for m in range(1,n+1):
  for pole in ['zero','infinity']:
   active=[j for j in range(k)if all(P[r,j]==0 for r in (range(m)if pole=='zero'else range(n-m+1,n+1)))];key=tuple(active)
   if not active or len(active)==k or key in groups:continue
   groups.append(key)
 out=[]
 for active in groups:
  h=len(active)
  if h==1:
   a=np.zeros(k,complex);a[active[0]]=1;b=npbary(a,V,P);out.append(dict(active=active,distance=float(np.linalg.norm(b)),basisCoordinates=[[z.real,z.imag]for z in a]));continue
  for chart in range(h):
   def unpack(x):
    a=np.zeros(k,complex);a[list(active)]=np.insert(x[:h-1]+1j*x[h-1:],chart,1);return a/np.linalg.norm(a)
   def obj(x):
    b=npbary(unpack(x),V,P);return -np.dot(b,b)
   res=differential_evolution(obj,[(-3,3)]*(2*h-2),seed=4700+n+h+chart,maxiter=600,popsize=15,tol=1e-10,polish=False)
   res2=minimize(obj,res.x,method='Nelder-Mead',options={'maxiter':5000,'fatol':1e-14,'xatol':1e-11});a=unpack(res2.x);b=npbary(a,V,P)
   out.append(dict(active=active,chart=chart,distance=float(np.linalg.norm(b)),basisCoordinates=[[z.real,z.imag]for z in a]));print(n,active,chart,'d',np.linalg.norm(b),flush=True)
 (ROOT/f'n{n}-poles.json').write_text(json.dumps(dict(best=max(out,key=lambda r:r['distance']),attempts=out),indent=2))
if __name__=='__main__':
 for n in map(int,sys.argv[1:]):run(n)
