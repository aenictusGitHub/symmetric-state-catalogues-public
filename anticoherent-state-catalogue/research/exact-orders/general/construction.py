"""Rational construction of pure states of exactly AC order t.

For every integer t>=1, positivity is guaranteed for N>=4*t*(t+1)**4.
The same formula is usable below that bound if every weight is positive.
Requires SymPy. No optimization or approximate linear algebra is used.
"""
from pathlib import Path
import argparse,json
import sympy as S

def power_sum(n,k):
    return S.Integer(n+1) if k==0 else S.cancel((S.bernoulli(k+1,n+1)-S.bernoulli(k+1,0))/(k+1))

def weight_polynomial(n,t):
    if t<1 or n<1:raise ValueError('N and t must be positive integers')
    s=t+1;m,remainder=divmod(n,s)
    if m<t:raise ValueError('This grid has fewer than t+1 nodes; another construction may exist')
    target=S.Matrix([power_sum(n,q)/(S.Integer(n+1)*n**q) for q in range(t+1)])
    moments=[S.cancel((s*S.Rational(s,n)**q*power_sum(m-1,q)+(remainder+1)*S.Rational(s*m,n)**q)/(n+1)) for q in range(2*t+1)]
    gram=S.Matrix(t+1,t+1,lambda a,b:moments[a+b])
    return [S.factor(x) for x in gram.inv()*target]

def weights(n,t):
    alpha=weight_polynomial(n,t);s=t+1;m,remainder=divmod(n,s)
    for i in range(m+1):
        x=S.Rational(s*i,n);q=S.Rational(s if i<m else remainder+1,n+1)
        w=S.factor(q*sum(a*x**j for j,a in enumerate(alpha)))
        if w<=0:raise ValueError('This grid has a nonpositive weight; theorem threshold is sufficient, not necessary')
        yield s*i,w

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('N',type=int);parser.add_argument('t',type=int);parser.add_argument('--output',type=Path)
    args=parser.parse_args();n,t=args.N,args.t
    support=list(weights(n,t))
    record={'n':n,'exactACOrder':t,'guaranteedThreshold':4*t*(t+1)**4,'weightPolynomialCoefficients':list(map(str,weight_polynomial(n,t))),'components':[{'r':r,'weight':str(w),'coefficient':str(S.sqrt(w))} for r,w in support],'phaseConvention':'All occupied Dicke coefficients are positive real.','nextOrderFailure':'The (0,t+1) entry of rho_(t+1) is a sum of strictly positive terms.'}
    text=json.dumps(record,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()
