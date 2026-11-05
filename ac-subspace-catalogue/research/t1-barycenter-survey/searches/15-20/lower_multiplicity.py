from search import *
for n in range(15,21):
 s=next(x for x in DATA['subspaces']if x['n']==n);k=s['dimension'];V=np.array([[complex(*s['basis'][j]['coefficients'][r])for j in range(k)]for r in range(n+1)]);P=V*np.array([(-1)**r*math.sqrt(math.comb(n,r))for r in range(n+1)])[:,None];m=k-2;rows=[]
 for chart in [0,1]:
  F=P if chart==0 else P[::-1];occupied=np.flatnonzero(np.linalg.norm(F,axis=1)>1e-15);lo,hi=occupied[0],occupied[-1]
  derivatives=[np.array([math.factorial(r)/math.factorial(r-l)*F[r]for r in range(l,n+1)])for l in range(m)]
  def evaluate(x):
   z=complex(*x[:2]);w=complex(*x[2:]);bad=(-1,np.zeros(k),np.zeros(3),z)
   if abs(z)<1e-4:return bad
   A=np.array([np.polynomial.polynomial.polyval(z,f)for f in derivatives]);A/=np.linalg.norm(A,axis=1)[:,None];Q=np.linalg.svd(A)[2][-2:].conj().T;c=Q@np.array([1,w]);c/=np.linalg.norm(c);q=(F@c)[lo:hi+1]
   for _ in range(m):
    q,r=np.polynomial.polynomial.polydiv(q,[-z,1])
    if np.linalg.norm(r)>1e-7:return bad
   rr=np.r_[[z]*m,np.roots(q[::-1])]
   if len(rr)+lo+n-hi!=n:return bad
   rad=abs(rr)**2;den=1+rad;b=np.array([sum(2*rr.real/den),sum(-2*rr.imag/den),sum((rad-1)/den)+n-hi-lo])/n
   if chart:b[1:]*=-1
   return np.linalg.norm(b),c,b,z
  def obj(x):return -evaluate(x)[0]
  for seed in range(2):
   opt=differential_evolution(obj,[(-1.2,1.2)]*2+[(-4,4)]*2,seed=9700+n*10+seed,popsize=16,maxiter=700,tol=1e-9,polish=True);d,c,b,z=evaluate(opt.x);rec={'n':n,'chart':chart,'seed':seed,'distance':float(d),'basisCoordinates':[[a.real,a.imag]for a in c],'barycenter':b.tolist(),'root':[z.real,z.imag],'rootMultiplicity':m,'evaluations':int(opt.nfev)};rows.append(rec);print(n,chart,seed,d,flush=True)
 (ROOT/f'n{n}-lower-multiplicity.json').write_text(json.dumps({'n':n,'attempts':4,'results':rows,'best':max(rows,key=lambda r:r['distance'])},indent=2)+'\n')
