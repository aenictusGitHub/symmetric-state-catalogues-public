#!/usr/bin/env python3
"""Exact higher-order witnesses in previously unresolved full partitions."""
import importlib.util
import json
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
def module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

verify=module('verify_states').verify
analyze=module('verify_multiplicities').analyze
definitions=[
    (8,{0:s.sqrt(s.Rational(7,27)),3:s.I*s.sqrt(s.Rational(4,27)),6:s.sqrt(s.Rational(16,27))},[2]+[1]*6,
     'Phase-shifted three-term state with one double point'),
    (10,{0:s.sqrt(s.Rational(1,6)),3:s.sqrt(s.Rational(1,3)),8:s.sqrt(s.Rational(1,2))},[2]+[1]*8,
     'Three-term state with one double point')]
states=[];checks=[]
for n,amplitudes,partition,label in definitions:
    c=[s.simplify(amplitudes.get(r,s.S.Zero)) for r in range(n+1)]
    state={'id':f'partition-n{n}-single-double-t2-sparse','n':n,'achievedOrder':2,
        'witnessesOrders':[1,2],'label':label,'status':'proven-existence',
        'support':sorted(amplitudes),'exactCoefficients':[s.sstr(v) for v in c],
        'coefficients':[[float(s.re(v)),float(s.im(v))] for v in c],
        'normalization':'Exact normalization and reduced matrices, independently recomputable.',
        'coefficientConvention':'Normalized Dicke basis; r counts excitations.',
        'references':[],'provenance':'Independent catalogue construction, 2026-10-28.',
        'notes':'Support gaps of at least three kill every off-diagonal condition through order two. The probabilities match the first two moments of the uniform Dicke distribution. Exact squarefree factorization proves one double point and all remaining points simple.'}
    exact=verify(state);roots=analyze(state)
    assert roots['partition']==partition
    checks.append({'id':state['id'],'partition':partition,'AC':exact,'multiplicity':roots,
        'probabilities':{str(r):s.sstr(s.simplify(c[r]*s.conjugate(c[r]))) for r in amplitudes},
        'previousLowerBound':1,'newLowerBound':2,'rigorousUpperBound':2})
    states.append(state)
output={'method':'Exact normalized-Dicke reduced matrices and algebraic squarefree factorization. No numerical optimization is used as proof.',
    'allVerifiedExactly':True,'states':states,'checks':checks}
(HERE/'ac-order-sparse.json').write_text(json.dumps(output,indent=2)+'\n')
print('Verified two exact improvements from 1–2 to maximum 2.')
