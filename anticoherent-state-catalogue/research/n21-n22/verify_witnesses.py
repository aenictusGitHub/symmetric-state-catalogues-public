"""Independent exact verification from the saved state components.

Uses a Schmidt matrix rather than the direct-entry sum in build_witnesses.py.
No numerical tolerance or restricted-state optimization is a proof premise.
Requires SymPy. Run from any directory.
"""
from pathlib import Path
import hashlib,json
import sympy as s

here=Path(__file__).resolve().parent
path=here/'states.json'
reports=[]
for state in json.loads(path.read_text())['states']:
    n,t=state['n'],state['achievedOrder']
    c=list(map(s.sympify,state['exactCoefficients']))
    assert len(c)==n+1
    assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
    reductions=[]
    for q in range(1,t+2):
        matrix=s.Matrix(q+1,n-q+1,lambda a,b:
            c[a+b]*s.sqrt(s.binomial(q,a)*s.binomial(n-q,b)/s.binomial(n,a+b)))
        delta=(matrix*matrix.H-s.eye(q+1)/(q+1)).applyfunc(s.simplify)
        zero=delta==s.zeros(q+1)
        assert zero==(q<=t)
        dist=s.simplify(s.trace(delta*delta.H))
        assert str(dist)==state['verification']['reductions'][q-1]['frobeniusDistanceSquared']
        reductions.append({'order':q,'equalsSymmetricIdentity':zero,'squaredFrobeniusDefect':str(dist)})
        if q==t+1:
            w=state['verification']['nextOrderWitness']
            assert s.simplify(delta[w['row'],w['column']]-s.sympify(w['differenceFromIdentity']))==0
    reports.append({'id':state['id'],'n':n,'exactOrder':t,'normalizationVerified':True,'reductions':reductions})
    print(f'PASS N={n}: exactly {t}-AC, not {t+1}-AC.')
(here/'independent-verification.json').write_text(json.dumps({'statesSha256':hashlib.sha256(path.read_bytes()).hexdigest(),'method':'Independent exact Schmidt matrices from the saved Dicke coefficients; every reduction through the next order.','records':reports,'passed':True},indent=2)+'\n')
