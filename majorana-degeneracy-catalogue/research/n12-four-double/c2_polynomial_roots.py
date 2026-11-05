#!/usr/bin/env python3
"""Direct seven-variable real-root search in independently derived C2 equations.

The polynomial equations are only necessary/sufficient within the stated C2
factor family. Every putative root is checked again with the normalized rho3.
No failure to find real roots is interpreted as a proof of nonexistence.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import root

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from search_high_order import FactorModel,separation,strings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--starts',type=int,default=1000)
    p.add_argument('--equations',type=Path,default=HERE/'c2-equations-input.json')
    p.add_argument('--output',type=Path,default=HERE/'c2-polynomial-root-results.json');a=p.parse_args()
    raw=a.equations.read_bytes();data=json.loads(raw);variables=sp.symbols(' '.join(data['variables']))
    expressions=[sp.sympify(s.replace('^','**'),locals=dict(zip(data['variables'],variables))) for s in data['equations']]
    scales=[max(abs(c) for c in sp.Poly(f,*variables).coeffs()) for f in expressions]
    F=sp.Matrix([f/s for f,s in zip(expressions,scales)])
    fun=sp.lambdify(variables,F,'numpy',cse=True);jac=sp.lambdify(variables,F.jacobian(variables),'numpy',cse=True)
    model=FactorModel([2]*4+[1]*4,3,{2:[0,2],1:[0,2]})
    assert np.linalg.norm(np.asarray(fun(1,-6,0,1,0,2,0),dtype=float))<1e-12
    rng=np.random.default_rng(12024413);runs=[];started=time.time()
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        for i in range(a.starts):
            x0=np.r_[np.exp(rng.normal(scale=.8)),rng.normal(scale=3,size=6)]
            if i%10==0:x0=np.array([1,-6,0,1,0,2,0])+rng.normal(scale=.1,size=7)
            solution=root(lambda x:np.asarray(fun(*x),dtype=float).ravel(),x0,
                          jac=lambda x:np.asarray(jac(*x),dtype=float),method='hybr',options={'maxfev':1200,'xtol':1e-11})
            x=solution.x
            if not np.all(np.isfinite(x)):
                runs.append({'index':i,'finite':False});continue
            aa,b,c,d,e,f,g=x
            factors=[np.array([aa,0,b+1j*c,0,1]),np.array([d+1j*e,0,f+1j*g,0,1])]
            model.evaluate(model.encode(factors));residual=float(np.linalg.norm(model.residual));gap,_=separation(factors)
            record={'index':i,'finite':True,'solverSuccess':bool(solution.success),'evaluations':solution.nfev,
                    'scaledPolynomialResidual':float(np.linalg.norm(np.asarray(fun(*x),dtype=float))),
                    'independentHSResidual':residual,'minimumDistinctRootSeparation':gap,
                    'maximumParameter':float(max(abs(x))),'realGaugeConstant':float(aa),
                    'provisionalSeparatedCandidate':residual<1e-10 and gap>.02}
            if record['provisionalSeparatedCandidate']:
                record['variables']=x.tolist();record['factors']=[strings(p) for p in factors]
                record['normalizedDickeCoefficients']=strings(model.coefficients)
            runs.append(record)
            if i%100==99 or i==a.starts-1:
                summary={'attempts':len(runs),'provisionalSeparatedCandidates':sum(r.get('provisionalSeparatedCandidate',False) for r in runs),
                         'smallResidualRoots':sum(r.get('independentHSResidual',1)<1e-10 for r in runs),
                         'bestHSResidualWithGapAbovePoint02':min((r['independentHSResidual'] for r in runs if r.get('minimumDistinctRootSeparation',0)>.02),default=None),
                         'runtimeSeconds':round(time.time()-started,3)}
                a.output.write_text(json.dumps({'scientificStatus':'Numerical real-root search only; no global exclusion or exact AC claim.',
                                               'equationInput':a.equations.name,'equationSHA256':hashlib.sha256(raw).hexdigest(),
                                               'seed':12024413,'summary':summary,'runs':runs},indent=2)+'\n')
                print(summary,flush=True)


if __name__=='__main__':main()
