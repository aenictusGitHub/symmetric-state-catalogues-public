"""Portable exact seven-anticoherent N=24 witness. Requires SymPy only."""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent

def schmidt_reduction(c,n,t):
    m=s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.Rational(
        s.binomial(t,a)*s.binomial(n-t,b),s.binomial(n,a+b))))
    return (m*m.H).applyfunc(s.simplify)

def main():
    n=24;z=s.symbols('z')
    f=s.Poly((1+14*z**4+z**8)**3,z)
    g=s.Poly(z**4*(1-z**4)**4,z)
    v=[f.nth(r)/s.sqrt(s.binomial(n,r)) for r in range(n+1)]
    w=[g.nth(r)/s.sqrt(s.binomial(n,r)) for r in range(n+1)]
    p,q,h=s.symbols('p q h',real=True)
    weights=[p*a*a+q*b*b+2*h*a*b for a,b in zip(v,w)]
    equations=[]
    for k in [0,2,4,6]:
        target=s.Rational(sum((2*r-n)**k for r in range(n+1)),n+1)
        equations.append(sum(ww*(2*r-n)**k for r,ww in enumerate(weights))-target)
    sol=s.solve(equations,(p,q,h))
    assert sol=={p:s.Rational(1131,16384),q:s.Rational(631449,320),h:s.Rational(22659,4096)}
    alpha=s.sqrt(sol[p])
    beta=s.simplify(sol[h]/alpha+s.I*s.sqrt(sol[q]-sol[h]**2/sol[p]))
    assert s.simplify(alpha-s.sqrt(1131)/128)==0
    assert s.simplify(beta-(581*s.sqrt(1131)/928+3289*s.sqrt(3045)*s.I/4640))==0
    c=[s.simplify(alpha*a+beta*b) for a,b in zip(v,w)]
    assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
    checks=[]
    for t in range(1,9):
        delta=(schmidt_reduction(c,n,t)-s.eye(t+1)/s.Integer(t+1)).applyfunc(s.simplify)
        zero=delta==s.zeros(t+1)
        assert zero==(t<=7)
        checks.append({'order':t,'equalsSymmetricIdentity':zero,
            'squaredFrobeniusDefect':str(s.simplify(s.trace(delta*delta.H)))})
        if t==8:assert delta[0,0]==s.Rational(13,155040)
    formula=(r'\alpha\bigl(|D_{24}^{(0)}\rangle+|D_{24}^{(24)}\rangle\bigr)'
        r'+\frac{42\alpha+\beta}{\sqrt{10626}}\bigl(|D_{24}^{(4)}\rangle+|D_{24}^{(20)}\rangle\bigr)'
        r'+\frac{591\alpha-4\beta}{\sqrt{735471}}\bigl(|D_{24}^{(8)}\rangle+|D_{24}^{(16)}\rangle\bigr)'
        r'+\frac{2828\alpha+6\beta}{\sqrt{2704156}}|D_{24}^{(12)}\rangle')
    state={'id':'octahedral-n24-t7','n':24,'achievedOrder':7,'witnessesOrders':list(range(1,8)),
        'label':'Exact octahedral state · order 7','status':'proven-existence',
        'stateLatex':formula,'parametersLatex':r'\alpha=\frac{\sqrt{1131}}{128},\quad\beta=\frac{581\sqrt{1131}}{928}+i\frac{3289\sqrt{3045}}{4640}',
        'support':[r for r,x in enumerate(c) if x!=0],
        'exactCoefficients':[str(x) for x in c],
        'coefficients':[[float(s.re(x)),float(s.im(x))] for x in c],
        'normalization':'sum_r |c_r|^2 = 1, verified exactly',
        'coefficientConvention':'Normalized Dicke basis; r counts excitations.',
        'provenance':'Exact octahedral binary forms with rational second-, fourth-, and sixth-moment matching.',
        'notes':'Exactly seventh-order anticoherent, not eighth-order. The global exclusion of eighth order is a separate bound.',
        'verification':{'method':'Exact Schmidt reductions through order eight.',
            'exactNormalization':True,'exactOrderVerified':True,'reductions':checks,
            'nextOrderWitness':{'order':8,'row':0,'column':0,'differenceFromIdentity':'13/155040'},
            'globalMaximumOrderProven':False}}
    (HERE/'states.json').write_text(json.dumps({'states':[state]},indent=2,ensure_ascii=False)+'\n')
    print('Built exact N=24 order-seven state; rho8[0,0]-1/9=13/155040.')

if __name__=='__main__':main()
