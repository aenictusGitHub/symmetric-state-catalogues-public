#!/usr/bin/env python3
"""Cyclic-support follow-up; numerical only, no exclusion from failure."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ.setdefault(key,'1')
from pathlib import Path
import json,time,itertools,argparse
import numpy as np
from scipy.optimize import least_squares
from search import problem,independent_check,HERE

def run(n,t,k,starts=20,steps=1500):
 rng=np.random.default_rng(4000000+n*1000+t*100+k);attempts=[];best=1.;bestB=None;t0=time.monotonic()
 patterns=[]
 for mod in range(max(k,2),t+3):
  for residues in itertools.combinations(range(mod),k):
   masks=np.column_stack([np.arange(n+1)%mod==r for r in residues]);patterns.append((mod,residues,masks))
 for mod,residues,mask in patterns:
  fun,jac,unpack,coords=problem(n,t,k,mask)
  for repeat in range(starts):
   B=(rng.normal(size=(n+1,k))+1j*rng.normal(size=(n+1,k)))*mask;B/=np.linalg.norm(B,axis=0)
   sol=least_squares(fun,coords(B),jac=jac,max_nfev=steps,ftol=1e-13,xtol=1e-13,gtol=1e-13)
   Q,check=independent_check(unpack(sol.x),n,t);res=check['maxCrossReductionFrobenius'];attempts.append({'modulus':mod,'residues':residues,'randomIndex':repeat,'nfev':sol.nfev,'residualL2':float(np.linalg.norm(sol.fun)),**check})
   if res<best:best=res;bestB=Q
   if res<1e-12:break
  print(n,t,k,mod,residues,best,round(time.monotonic()-t0,1),flush=True)
  if best<1e-12:break
 report={'n':n,'t':t,'dimension':k,'method':'Disjoint cyclic-support complex least squares with exact analytic Jacobian; independent physical reduced-matrix check.','seed':4000000+n*1000+t*100+k,'attempts':attempts,'elapsedSeconds':time.monotonic()-t0,'maxCrossReductionFrobenius':best,'found':best<1e-10,'bestBasis':[[[float(z.real),float(z.imag)]for z in bestB[:,i]]for i in range(k)],'scope':'Restricted numerical search. Failure excludes neither unrestricted subspaces nor all symmetric support patterns.'}
 (HERE/f'structured-n{n}-t{t}-k{k}.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--t',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--starts',type=int,default=20);p.add_argument('--steps',type=int,default=1500);run(**vars(p.parse_args()))
