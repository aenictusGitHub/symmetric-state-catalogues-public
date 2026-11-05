"""Interpolate a five-band cyclic solution to larger dimensions; exploratory only."""
from pathlib import Path
import json,argparse,time
import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline
from numeric_cyclic import problem
HERE=Path(__file__).resolve().parent

def run(N,k,tries=40):
 source=json.loads((HERE/'dihedral-n24-k5-sign1.json').read_text());z=np.array([complex(*v) for v in source['components']]);M=np.array([z[a::5] for a in range(5)])
 M*=np.exp(-1j*np.angle(M[:,0]))[:,None]
 xold=np.arange(5)/4;xnew=np.arange(k)/(k-1)
 amplitudes=CubicSpline(xold,np.log(abs(M)),axis=0)(xnew)
 phases=CubicSpline(xold,np.unwrap(np.angle(M),axis=0),axis=0)(xnew)
 C=np.exp(amplitudes+1j*phases);C/=np.linalg.norm(C,axis=1)[:,None]
 init=np.zeros(N+1,dtype=complex)
 for a in range(k):init[a::k]=C[a]
 fun,_,unpack=problem(N,k);rng=np.random.default_rng(449);best=1e3;t=time.time()
 for restart in range(tries):
  cc=init.copy()
  if restart:
   scale=.02*1.4**(restart%12);cc*=np.exp(rng.normal(size=N+1)*scale+1j*rng.normal(size=N+1)*scale)
  for a in range(k):cc[a::k]*=np.exp(-1j*np.angle(cc[a]));cc[a::k]/=np.linalg.norm(cc[a::k])
  x=np.r_[cc.real,cc[k:].imag]
  sol=least_squares(fun,x,jac=lambda x:fun(x,True),ftol=1e-13,gtol=1e-13,xtol=1e-13,max_nfev=2500)
  val=np.linalg.norm(sol.fun);print(N,k,restart,val,sol.nfev,round(time.time()-t,1),flush=True)
  if val<best:
   best=val;c=unpack(sol.x);(HERE/f'continuation-n{N}-k{k}.json').write_text(json.dumps({'N':N,'k':k,'residualNorm':float(val),'components':[[float(v.real),float(v.imag)] for v in c],'restart':restart,'source':'dihedral-n24-k5-sign1.json','note':'Numerical interpolation search; no exclusion proof.'},indent=2)+'\n')
  if best<1e-10:break
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--N',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--tries',type=int,default=40);a=p.parse_args();run(a.N,a.k,a.tries)
