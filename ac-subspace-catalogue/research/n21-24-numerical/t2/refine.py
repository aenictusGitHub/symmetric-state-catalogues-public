#!/usr/bin/env python3
"""Newton refinement of accepted floating-point bases; still numerical evidence."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json
import numpy as np
from search import problem,coords,check,basis_json
HERE=Path(__file__).resolve().parent
for n in [21,22,23]:
 path=HERE/f'n{n}-k4.json';d=json.loads(path.read_text());assert d['found'];B=np.array([[complex(*z)for z in b['coefficients']]for b in d['basis']]).T
 fun,jac,unpack,_=problem(n,4);x=coords(B);old=d['bestCrossReductionFrobenius'];history=[]
 for step in range(5):
  f=fun(x);J=jac(x);delta=np.linalg.lstsq(J,-f,rcond=1e-11)[0];x+=delta
  history.append({'constraintResidualL2':float(np.linalg.norm(fun(x))),'newtonStepNorm':float(np.linalg.norm(delta))})
 B=np.linalg.qr(unpack(x))[0];c=check(B,n);assert c['maxOrderedCrossReductionFrobenius']<old
 d['basis']=basis_json(B);d['refinement']={'method':'Five minimum-norm Newton steps in full complex Dicke coordinates, followed by QR.','previousMaxCrossReductionFrobenius':old,'steps':history};d['independentPhysicalCheck']=c;d['bestCrossReductionFrobenius']=c['maxOrderedCrossReductionFrobenius'];path.write_text(json.dumps(d,indent=2)+'\n');print(n,c['maxOrderedCrossReductionFrobenius'])
