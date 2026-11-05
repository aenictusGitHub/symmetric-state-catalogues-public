import json,math,numpy as np
from pathlib import Path

def pauli_sector(c,k):
 n=len(c)-1;idx=np.arange(2**n)
 pop=np.array([int(i).bit_count() for i in idx])
 psi=np.asarray([c[w]/math.sqrt(math.comb(n,w)) for w in pop]);total=0.
 for x in range(k+1):
  for y in range(k-x+1):
   z=k-x-y;flip=(1<<(x+y))-1
   yz=((1<<y)-1)<<x | ((1<<z)-1)<<(x+y)
   phase=(1j)**y*(-1.)**np.array([int(i&yz).bit_count() for i in idx])
   expect=np.vdot(psi[idx^flip],phase*psi)
   total+=math.comb(n,k)*math.comb(k,x)*math.comb(k-x,y)*expect.real**2
 return total
R=Path(__file__).resolve().parent;v=json.loads((R/'verification.json').read_text())
c=np.zeros(12,complex)
for r,x in v['coefficients100digits'].items():c[int(r)]=float(x)
t=np.zeros(12,complex)
for r,x in enumerate(c):t[11-r]=(-1)**r*x.conjugate()
assert abs(np.vdot(c,t))<1e-14
baseline=np.array([pauli_sector(c,k) for k in (2,4,6,8,10)])
rng=np.random.default_rng(2711);err=0.
for i in range(12):
 z=rng.normal(size=2)+1j*rng.normal(size=2);z/=np.linalg.norm(z)
 family=z[0]*c+z[1]*t
 assert abs(np.linalg.norm(family)-1)<1e-14
 sectors=np.array([pauli_sector(family,k) for k in (2,4,6,8,10)])
 err=max(err,float(np.max(abs(sectors-baseline))))
assert err<1e-9
v['kramersFamilyNumericalChecks']=12;v['kramersFamilyMaximumSectorError']=err
(R/'verification.json').write_text(json.dumps(v,indent=2)+'\n')
print('12 random complex Kramers mixtures preserve all five even sectors; maxerror=',err)
