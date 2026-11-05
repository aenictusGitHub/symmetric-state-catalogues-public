import sympy as S
import json
from pathlib import Path

def sums(n,q):
    if q==0:return S.Integer(n+1)
    return S.cancel((S.bernoulli(q+1,n+1)-S.bernoulli(q+1,0))/(q+1))

def construction(n,t):
    step=t+1;m,r=divmod(n,step)
    if m<t:raise ValueError('Too few grid points')
    mu=S.Matrix([sums(n,q)/(S.Integer(n+1)*n**q) for q in range(t+1)])
    moments=[S.cancel((step*(S.Rational(step,n))**q*sums(m-1,q)+(r+1)*S.Rational(step*m,n)**q)/(n+1)) for q in range(2*t+1)]
    G=S.Matrix(t+1,t+1,lambda a,b:moments[a+b]);alpha=G.inv()*mu
    return [S.factor(a) for a in alpha],moments,mu

def full_reduction(n,t):
    a,_,_=construction(n,t);step=t+1;m,r=divmod(n,step)
    c={step*i:S.sqrt(S.Rational(step if i<m else r+1,n+1)*sum(a[q]*S.Rational(step*i,n)**q for q in range(t+1))) for i in range(m+1)}
    assert all(S.simplify(x*x)>0 for x in c.values())
    out=[]
    for q in [t,t+1]:
        mat=S.Matrix(q+1,n-q+1,lambda i,j:c.get(i+j,0)*S.sqrt(S.Rational(S.binomial(q,i)*S.binomial(n-q,j),S.binomial(n,i+j))))
        delta=(mat*mat.T-S.eye(q+1)/(q+1)).applyfunc(S.simplify)
        assert all(x==0 for x in delta) if q==t else delta[0,q]>0
        out.append({'q':q,'maximallyMixed':q==t,'endpointCoherence':str(delta[0,q])})
    return {'n':n,'t':t,'reductions':out}

report={'smallExactChecks':[],'thresholdChecks':[]}
for n,t in [(2,1),(3,1),(6,2),(8,2),(12,3),(20,4),(30,5),(42,6)]:
    report['smallExactChecks'].append(full_reduction(n,t));print('Full reductions',n,t,'pass',flush=True)
for t in range(1,9):
    n=4*t*(t+1)**4;a,g,mu=construction(n,t)
    # Bernstein form is a distinct exact positivity check on the whole interval.
    beta=[S.factor(sum(a[j]*S.binomial(k,j)/S.binomial(t,j) for j in range(k+1))) for k in range(t+1)]
    assert min(beta)>0
    assert all(S.cancel(sum(a[j]*g[q+j] for j in range(t+1))-mu[q])==0 for q in range(t+1))
    report['thresholdChecks'].append({'n':n,'t':t,'allBernsteinCoefficientsPositive':True,'minimumBernsteinCoefficient':str(min(beta))})
    print('Threshold',n,t,'positive',flush=True)
# Independent symbolic verification of the compact t=2 family.
n,k=S.symbols('N k',positive=True)
w=[(3*k-n+1)/(6*k), n*(n-1)/(6*k*(n-k)),(2*n+1-3*k)/(6*(n-k))]
assert S.simplify(sum(w)-1)==0
assert S.simplify(w[1]*k+w[2]*n-n/2)==0
assert S.simplify(w[1]*k*k+w[2]*n*n-n*(2*n+1)/6)==0
assert S.simplify(w[1]*k**3+w[2]*n**3-n*n*(n+1)/4-n*(n-1)*(-2*k+n)/12)==0
report['threeTermOrderTwoSymbolicIdentities']=True
Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')
print('ALL CHECKS PASS')
