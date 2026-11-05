#!/usr/bin/env python3
"""Long full-complex continuation of the best unsuccessful dimension-five fit."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
from scipy.optimize import least_squares
from search import problem,coords,basis_json,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,choices=[22,23],required=True);p.add_argument('--steps',type=int,default=5000);a=p.parse_args();path=HERE/f'n{a.n}-k5.json';raw=path.read_bytes();d=json.loads(raw);B=np.array([[complex(*z)for z in r['coefficients']]for r in d['basis']]).T
 f,j,u,_=problem(a.n,5);start=time.monotonic();sol=least_squares(f,coords(B),jac=j,max_nfev=a.steps,ftol=1e-13,gtol=1e-13,xtol=1e-13)
 C=np.linalg.qr(u(sol.x))[0];physical=check(C,a.n);value=physical['maxOrderedCrossReductionFrobenius']
 item={'seedName':'long-continuation-of-best-full-complex-fit','inputRecordSHA256':hashlib.sha256(raw).hexdigest(),'nfev':int(sol.nfev),'optimizerStatus':int(sol.status),'constraintResidualL2':float(np.linalg.norm(sol.fun)),'maxCrossReductionFrobenius':value,'orthonormalityFrobenius':physical['orthonormalityFrobenius'],'elapsedSeconds':time.monotonic()-start}
 d['attempts'].append(item);d['totalElapsedSeconds']+=item['elapsedSeconds']
 if value<d['bestCrossReductionFrobenius']:d['basis']=basis_json(C);d['bestCrossReductionFrobenius']=value;d['independentPhysicalCheck']=physical;d['bestSeed']=item['seedName']
 d['found']=d['bestCrossReductionFrobenius']<1e-10;path.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(item),flush=True)
