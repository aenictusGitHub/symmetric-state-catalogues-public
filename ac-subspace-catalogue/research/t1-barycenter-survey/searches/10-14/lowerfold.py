from survey import *
def run(n):
 s,V,M,Pmp=get(n);P=np.array([[complex(x)for x in row]for row in Pmp]);k=s['dimension'];m=k-1
 def getrec(root):
  D=np.array([sum(P[r]*math.factorial(r)/math.factorial(r-l)*root**(r-l)for r in range(l,n+1))for l in range(m)]);D/=np.linalg.norm(D,axis=1)[:,None]
  _,sv,vh=np.linalg.svd(D);a=vh[-1].conj();p=P@a;hi=max(r for r in range(n+1)if np.any(P[r]!=0));p=p[:hi+1]
  q,rem=np.polynomial.polynomial.polydiv(p,np.polynomial.polynomial.polypow([-root,1],m));rr=np.r_[np.repeat(root,m),np.roots(q[::-1])];rad=abs(rr)**2;b=np.array([sum(2*rr.real/(1+rad)),sum(-2*rr.imag/(1+rad)),sum((rad-1)/(1+rad))+n-hi])/n
  return b,a,np.linalg.norm(rem)/np.linalg.norm(p)
 def objective(x):
  b,a,err=getrec(complex(*x))
  return -np.dot(b,b)if err<1e-7 else 0
 out=[]
 for seed in range(8):
  res=differential_evolution(objective,[(-4,4)]*2,seed=7000+n*10+seed,popsize=20,maxiter=450,tol=1e-11,polish=True);b,a,err=getrec(complex(*res.x))
  row=dict(root=res.x.tolist(),distance=float(np.linalg.norm(b)),basisCoordinates=[[z.real,z.imag]for z in a],remainder=err);out.append(row);print(n,seed,row,flush=True)
 (ROOT/f'n{n}-lowerfold.json').write_text(json.dumps(dict(multiplicity=m,attemptCount=len(out),attempts=out,best=max(out,key=lambda r:r['distance'])),indent=2))
if __name__=='__main__':
 for n in map(int,sys.argv[1:]):run(n)
