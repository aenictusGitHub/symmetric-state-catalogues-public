"""Release cyclic support constraints in a numerical saturation search."""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import least_squares
from general_search import problem
HERE=Path(__file__).resolve().parent
n,k=29,6;fun,jac,unpack=problem(n,k)
d=json.loads((HERE/'cyclic-n29-k6.json').read_text());z=np.array([complex(*c) for c in d['components']]);B=np.zeros((n+1,k),complex)
for a in range(k):B[a::k,a]=z[a::k]
rng=np.random.default_rng(2410);best=1e9;report=[];start=time.time()
for noise in [0.001,0.03]:
 C=B+noise*(rng.normal(size=B.shape)+1j*rng.normal(size=B.shape));x=np.r_[C.real.reshape(-1),C.imag.reshape(-1)]
 sol=least_squares(fun,x,jac=jac,max_nfev=2500,ftol=1e-13,gtol=1e-13,xtol=1e-13)
 norm=float(np.linalg.norm(sol.fun));print(noise,sol.nfev,norm,time.time()-start,flush=True)
 report.append({'perturbation':noise,'evaluations':sol.nfev,'residual':norm})
 if norm<best:
  best=norm;C=unpack(sol.x)
  (HERE/'unrestricted-from-cyclic-n29-k6.json').write_text(json.dumps({'n':n,'k':k,'residual':norm,'basis':[[[float(z.real),float(z.imag)] for z in C[:,a]] for a in range(k)],'trials':report,'note':'Numerical result only; neither existence nor nonexistence certified.'},indent=2)+'\n')
