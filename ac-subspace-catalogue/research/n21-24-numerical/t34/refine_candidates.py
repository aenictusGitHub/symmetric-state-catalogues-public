#!/usr/bin/env python3
"""Refine saved candidates and check them through physical reduced matrices."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ.setdefault(key,'1')
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from search import problem,independent_check,HERE

def main():
 records=[]
 for file in sorted(HERE.glob('n*-t*-k*.json')):
  d=json.loads(file.read_text())
  if not d['found']:continue
  n,t,k=d['n'],d['t'],d['dimension'];B=np.array([[complex(*z)for z in row]for row in d['bestBasis']]).T
  fun,jac,unpack,coords=problem(n,t,k)
  sol=least_squares(fun,coords(B),jac=jac,max_nfev=1500,ftol=1e-15,xtol=1e-15,gtol=1e-15)
  Q,check=independent_check(unpack(sol.x),n,t);assert check['maxCrossReductionFrobenius']<1e-12
  rec={'id':f'numerical-search-n{n}-k{k}-t{t}','n':n,'t':t,'dimension':k,'label':f'Numerically optimized AC{t} subspace on {n} qubits','basis':[{'coefficients':[[float(z.real),float(z.imag)]for z in Q[:,i]]}for i in range(k)],'verification':{'method':'Independent normalized-Dicke splitting after QR; all ordered cross-reduction matrices checked.','nfevRefinement':int(sol.nfev),'refinedTensorResidualL2':float(np.linalg.norm(sol.fun)),**check},'searchSource':file.name,'evidence':'numerical','notes':'High-accuracy numerical whole-subspace candidate. This does not establish exact existence or global optimality.'}
  records.append(rec);print(n,t,k,check['maxCrossReductionFrobenius'],flush=True)
 (HERE/'candidates.json').write_text(json.dumps({'candidates':records},indent=2)+'\n')
if __name__=='__main__':main()
