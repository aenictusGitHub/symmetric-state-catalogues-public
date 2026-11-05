"""Independently verify saved exact-order states by a Schmidt matrix.

No use of the construction search or its reduced-state implementation.
Requires SymPy. Run from any working directory.
"""
from pathlib import Path
import json
import sympy as S

HERE=Path(__file__).resolve().parent
states=json.loads((HERE/'states.json').read_text())['states']
checks=[]
for state in states:
    n,t=state['n'],state['achievedOrder']
    c=list(map(S.sympify,state['exactCoefficients']))
    assert S.simplify(sum(x*S.conjugate(x) for x in c))==1
    for q in range(1,t+2):
        # |D_n^r> = sum_a sqrt(C(q,a)C(n-q,r-a)/C(n,r)) |D_q^a>|D_(n-q)^(r-a)>.
        schmidt=S.Matrix(q+1,n-q+1,lambda a,b:c[a+b]*S.sqrt(S.Rational(S.binomial(q,a)*S.binomial(n-q,b),S.binomial(n,a+b))))
        delta=(schmidt*schmidt.conjugate().T-S.eye(q+1)/(q+1)).applyfunc(S.simplify)
        zero=all(x==0 for x in delta)
        assert zero==(q<=t),(state['id'],q)
        stored=next(r for r in state['verification']['reductions'] if r['t']==q)
        assert S.simplify(sum(S.expand_complex(x*S.conjugate(x)) for x in delta)-S.sympify(stored['frobeniusDistanceSquared']))==0
        if q==t+1:
            witness=state['verification']['nextOrderWitness']
            assert delta[witness['row'],witness['column']]!=0
            assert S.simplify(delta[witness['row'],witness['column']]-S.sympify(witness['differenceFromIdentity']))==0
    checks.append({'id':state['id'],'n':n,'exactOrder':t,'passed':True})
result={'method':'Independent Schmidt-matrix partial trace, exact SymPy arithmetic','passed':True,'statesVerified':len(checks),'checks':checks}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'PASS: {len(checks)} normalized exact-order witnesses; all lower reductions and every next-order obstruction verified.')
