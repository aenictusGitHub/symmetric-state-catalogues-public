"""Independent direct partial trace from saved coefficients. SymPy only."""
from pathlib import Path
import hashlib,json
import sympy as s

def main():
    here=Path(__file__).resolve().parent;path=here/'states.json'
    states=json.loads(path.read_text())['states'];reports=[]
    for state in states:
        n,t=state['n'],state['achievedOrder'];c=list(map(s.sympify,state['exactCoefficients']))
        assert n==24 and t==7 and len(c)==25
        assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
        support=[r for r,x in enumerate(c) if x!=0]
        checks=[]
        for q in range(1,t+2):
            rho=s.zeros(q+1)
            for a in range(q+1):
                for b in range(a,q+1):
                    value=0
                    for r in support:
                        ell=r-a
                        if 0<=ell<=n-q and 0<=b+ell<=n:
                            value+=c[r]*s.conjugate(c[b+ell])*s.binomial(n-q,ell)*s.sqrt(s.Rational(
                                s.binomial(q,a)*s.binomial(q,b),s.binomial(n,r)*s.binomial(n,b+ell)))
                    rho[a,b]=s.simplify(value);rho[b,a]=s.conjugate(rho[a,b])
            delta=(rho-s.eye(q+1)/s.Integer(q+1)).applyfunc(s.simplify)
            zero=delta==s.zeros(q+1);assert zero==(q<=t)
            dist=s.simplify(s.trace(delta*delta.H))
            assert s.simplify(dist-s.sympify(state['verification']['reductions'][q-1]['squaredFrobeniusDefect']))==0
            checks.append({'order':q,'equalsSymmetricIdentity':zero,'squaredFrobeniusDefect':str(dist)})
            if q==8:assert delta[0,0]==s.Rational(13,155040)
        z=s.symbols('z')
        majorana=s.Poly(sum((-1)**r*s.sqrt(s.binomial(n,r))*x*z**r
            for r,x in enumerate(c)),z,extension=True)
        assert majorana.degree()==n and c[0]!=0
        assert s.gcd(majorana,majorana.diff()).degree()==0
        reports.append({'id':state['id'],'n':n,'exactOrder':t,'normalizationVerified':True,
            'majoranaStarsAllDistinct':True,'reductions':checks})
    report={'passed':True,'statesSha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'method':'Independent direct Dicke partial-trace sums from saved coefficients; no construction imports.',
        'records':reports}
    (here/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: exactly seventh-order N=24 state, all reductions independently verified.')

if __name__=='__main__':main()
