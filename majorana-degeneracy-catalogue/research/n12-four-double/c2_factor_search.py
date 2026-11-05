#!/usr/bin/env python3
"""Independent starts in the complete complex even-quartic Q^2 R family.

This is a C2 symmetry restriction, not a proof that every possible state has C2.
The reused factor model computes normalized states and analytic real derivatives.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from search_high_order import FactorModel,separation,strings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--starts',type=int,default=1000)
    p.add_argument('--output',type=Path,default=HERE/'c2-factor-results.json')
    a=p.parse_args();rng=np.random.default_rng(120442)
    model=FactorModel([2]*4+[1]*4,3,{2:[0,2],1:[0,2]})
    records=[];started=time.time()
    for index in range(a.starts):
        factors=[]
        for _ in range(2):
            roots=np.exp(rng.normal(scale=.6,size=2)+1j*rng.uniform(0,2*np.pi,size=2))
            factors.append(np.array([roots[0]*roots[1],0,-sum(roots),0,1]))
        x=model.encode(factors)
        result=least_squares(model.fun,x,jac=model.jac,max_nfev=1200,ftol=2e-14,xtol=2e-14,gtol=2e-14)
        model.evaluate(result.x);factors=model.factors(result.x)
        gap,_=separation(factors);residual=float(np.linalg.norm(model.residual));c=model.coefficients
        record={'index':index,'HSResidual':residual,'minimumDistinctRootSeparation':gap,
                'provisionalSeparatedCandidate':residual<1e-11 and gap>.02,
                'timeReversalInfidelity':float(max(0,1-abs(np.vdot(c,(-1)**np.arange(13)*np.conj(c[::-1])))**2)),
                'evaluations':result.nfev}
        if record['provisionalSeparatedCandidate']:
            record['factors']=[strings(f) for f in factors]
            record['normalizedDickeCoefficients']=strings(c)
        records.append(record)
        if index%50==49 or index==a.starts-1:
            summary={'attempts':len(records),'provisionalSeparatedCandidates':sum(r['provisionalSeparatedCandidate'] for r in records),
                     'smallResidualWithGapAbovePoint02':min((r['HSResidual'] for r in records if r['minimumDistinctRootSeparation']>.02),default=None),
                     'minimumHSResidual':min(r['HSResidual'] for r in records),'runtimeSeconds':round(time.time()-started,3)}
            a.output.write_text(json.dumps({'scientificStatus':'Numerical C2-restricted search. No exact existence or global exclusion is inferred.',
                                           'seed':120442,'family':'(z4+a*z2+b)^2*(z4+c*z2+d), a,b,c,d complex',
                                           'summary':summary,'runs':records},indent=2)+'\n')
            print(summary,flush=True)


if __name__=='__main__':main()
