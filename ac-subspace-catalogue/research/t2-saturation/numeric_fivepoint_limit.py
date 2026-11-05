"""Check the semiclassical five-band weighted spherical 2-design ansatz numerically."""
import numpy as np
from scipy.optimize import least_squares

def run(x,restarts=30):
 z=1-2*(np.arange(5)+x)/5;h=np.sqrt(1-z*z)
 def fun(y):
  w=y[:5];phi=np.r_[0,y[5:]];a=np.exp(1j*phi)
  b=np.array([sum(w*h*a),sum(w*z*h*a),sum(w*h*h*a*a)])
  return np.r_[sum(w)-1,sum(w*z),sum(w*z*z)-1/3,b.real,b.imag]
 rng=np.random.default_rng(657);best=10
 for i in range(restarts):
  sol=least_squares(fun,np.r_[np.ones(5)/5,rng.uniform(-np.pi,np.pi,4)],bounds=(np.r_[np.ones(5)*1e-8,np.ones(4)*-np.inf],np.inf),gtol=1e-12,xtol=1e-12,ftol=1e-12,max_nfev=2000)
  best=min(best,np.linalg.norm(sol.fun))
  if best<1e-10:break
 print(x,best,sol.x[:5],flush=True)
for x in [.01,.1,.25,.4,.5,.6,.75,.9,.99]:run(x)
