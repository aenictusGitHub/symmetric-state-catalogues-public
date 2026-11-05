#!/usr/bin/env python3
"""C2 searches including a pair of fixed pole stars (simple or double)."""
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from search_high_order import FactorModel,strings


class PoleModel(FactorModel):
    def __init__(self,pole_group):
        self.pole_group=pole_group
        allowed={2:[1],1:[0,2]} if pole_group==0 else {2:[0,2],1:[1]}
        super().__init__([2]*4+[1]*4,3,allowed)
    def factors(self,x):
        factors=super().factors(x)
        factors[self.pole_group][4]=0
        factors[self.pole_group][3]=1
        return factors


def root_gap(factors):
    roots=[z for p in factors for z in np.roots(np.trim_zeros(p[::-1],'f'))]
    # Exactly one factor has one projective root at infinity.
    distances=[2/math.sqrt(1+abs(z)**2) for z in roots]
    distances.extend(2*abs(a-b)/math.sqrt((1+abs(a)**2)*(1+abs(b)**2))
                     for i,a in enumerate(roots) for b in roots[:i])
    return float(min(distances))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--starts',type=int,default=300)
    p.add_argument('--output',type=Path,default=HERE/'c2-pole-results.json');a=p.parse_args()
    rng=np.random.default_rng(12024412);records=[]
    for pole_group in (0,1):
        model=PoleModel(pole_group)
        for i in range(a.starts):
            values=rng.normal(size=3)+1j*rng.normal(size=3)
            x=np.r_[values.real,values.imag]
            result=least_squares(model.fun,x,jac=model.jac,max_nfev=1500,ftol=2e-14,xtol=2e-14,gtol=2e-14)
            model.evaluate(result.x);factors=model.factors(result.x);gap=root_gap(factors)
            residual=float(np.linalg.norm(model.residual))
            record={'branch':'double-poles' if pole_group==0 else 'simple-poles','index':i,
                    'HSResidual':residual,'minimumProjectiveRootSeparation':gap,
                    'provisionalSeparatedCandidate':residual<1e-11 and gap>.02,'evaluations':result.nfev}
            if record['provisionalSeparatedCandidate']:
                record['factors']=[strings(f) for f in factors]
                record['normalizedDickeCoefficients']=strings(model.coefficients)
            records.append(record)
        print(pole_group,'best',min(r['HSResidual'] for r in records if r['branch']==record['branch']),
              'candidates',sum(r['provisionalSeparatedCandidate'] for r in records),flush=True)
        a.output.write_text(json.dumps({'scientificStatus':'Numerical fixed-pole C2 searches only; no existence or nonexistence proof.',
                                       'seed':12024412,'summary':{'attempts':len(records),'provisionalSeparatedCandidates':sum(r['provisionalSeparatedCandidate'] for r in records)},
                                       'runs':records},indent=2)+'\n')


if __name__=='__main__':main()
