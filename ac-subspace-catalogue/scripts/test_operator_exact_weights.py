#!/usr/bin/env python3
"""Exact arithmetic and independent Clebsch-Gordan regression checks.
Run: python3 scripts/test_operator_exact_weights.py (requires SymPy).
"""
from pathlib import Path
import hashlib,json
import sympy as s
from sympy.physics.wigner import clebsch_gordan
from build_operator_exact_weights import parse_exact,tensor,state_weights,Unsupported
ROOT=Path(__file__).resolve().parents[1]
# Arbitrary Python and decimal approximations are never accepted as exact input.
for bad in ['0.1','__import__("os")','unknown_parameter','sqrt(2,3)']:
    try:parse_exact(bad)
    except (Unsupported,SyntaxError):pass
    else:raise AssertionError(bad)
assert parse_exact('sqrt(2)/2 + I/3')==s.sqrt(2)/2+s.I/3
# Independent factorial-CG tensor coefficients fix phase, index and normalization.
for n,l,m,r in [(1,1,0,0),(4,2,1,1),(7,6,5,0),(8,4,4,2),(13,8,4,5),(26,2,0,9),(42,6,3,7)]:
    j=s.Rational(n,2);c=r+m
    expected=(-1)**c*clebsch_gordan(j,j,l,j-r,-j+c,m)
    assert s.simplify(tensor(n,l,m,r)-expected)==0,(n,l,m,r)
records={r['id']:r for r in json.loads((ROOT/'data/catalogue.json').read_text())['subspaces']}
fixtures={'radical-n13-k2-t3':{8:s.Rational(1365,13376)},'acliege-Π_7_2_2':{6:s.Rational(3,8)},'line-ac-n8-t3':{4:s.Rational(98,429)},'reflection-pairs-n26-k9-t1':{2:s.Rational(242,356265)}}
for id,weights in fixtures.items():
    result=state_weights(records[id])
    for rank,square in weights.items():assert parse_exact(result['weights'][str(rank)]['squared'])==square
payload=json.loads((ROOT/'data/operator-exact-weights.json').read_text())
assert payload['sourceSha256']==hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest()
for id,state in payload['states'].items():
    assert records[id]['evidence']=='exact'
    assert state['basisOrthonormality']=='exact' and state['multipolePurity']=='exact'
    for rank,w in state['weights'].items():
        square=parse_exact(w['squared'])
        assert abs(float(s.N(s.sqrt(square),17))-w['numericValue'])<1e-15
        assert len(w['latex'])<=80
print(json.dumps({'exactStates':len(payload['states']),'weights':sum(len(r['weights']) for r in payload['states'].values()),'independentCG':'PASS','sourceHash':'PASS','decimalAndCodeRejection':'PASS','exactFixtures':'PASS'}))
