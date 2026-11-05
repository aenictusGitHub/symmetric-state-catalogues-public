#!/usr/bin/env python3
"""Search antipodal axes with multiplicities(2,2,1,1); numerical results only."""
from pathlib import Path
import json,time,math
import numpy as np
from scipy.optimize import least_squares
N=12
bc=np.sqrt([math.comb(N,r) for r in range(N+1)])
idx=np.array([[a+b for b in range(11)] for a in range(3)])
w=np.array([[math.sqrt(math.comb(2,a)*math.comb(10,b)/math.comb(12,a+b)) for b in range(11)] for a in range(3)])
def q(th,ph):return np.array([-.5*np.sin(th)*np.exp(-1j*ph),np.cos(th),.5*np.sin(th)*np.exp(1j*ph)])
def coeff(x):
 p=np.array([0,0,1],complex)
 for f in [q(x[0],0),q(x[0],0),q(x[1],x[2]),q(x[3],x[4])]:p=np.convolve(p,f)
 c=np.pad(p,(0,13-len(p)))/bc
 return c/np.linalg.norm(c)
def residual(x):
 c=coeff(x);A=c[idx]*w;r=A@A.conj().T-np.eye(3)/3
 return np.array([r[0,0].real,r[1,1].real,r[2,2].real,*[y for a,b in [(1,0),(2,0),(2,1)] for y in [2**.5*r[a,b].real,2**.5*r[a,b].imag]]])
def axes(x):
 return [[0.,0.,1.]]+[[np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)] for t,p in [(x[0],0),(x[1],x[2]),(x[3],x[4])]]
def gap(x):
 a=axes(x);return min(np.sqrt(max(0.,2-2*abs(np.dot(v,u)))) for i,v in enumerate(a) for u in a[:i])
def main():
 rng=np.random.default_rng(124413);out=[];start=time.time()
 for i in range(1500):
  x=rng.uniform(-np.pi,np.pi,5)
  z=least_squares(residual,x,xtol=2e-13,ftol=2e-13,gtol=2e-13,max_nfev=350)
  rr=np.linalg.norm(z.fun);g=gap(z.x)
  if i<10 or (rr<1e-8 and g>.03) or i%100==0:print(i,'r',rr,'gap',g,flush=True)
  out.append({'start':i,'residual':rr,'minAxisSeparation':g,'angles':z.x.tolist(),'axes':axes(z.x),'coefficients':[[v.real,v.imag] for v in coeff(z.x)]})
  if rr<1e-10 and g>.05:print('CANDIDATE',out[-1],flush=True);break
 result={'method':'Finite multistart numerical search in general antipodal four-axis family. Not an exclusion or exact existence proof.','runs':out,'elapsed':time.time()-start}
 Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
