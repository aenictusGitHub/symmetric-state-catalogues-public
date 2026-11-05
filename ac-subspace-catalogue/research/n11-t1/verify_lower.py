#!/usr/bin/env python3
"""Check the displayed N11 dimension-four witness using exact Dicke splitting."""
from pathlib import Path
import json
import sympy as s

def verify():
    alpha=(36+5*s.I*s.sqrt(789))/(70*s.sqrt(11))
    supports=[{4:s.sqrt(s.Rational(7,10)),9:s.sqrt(s.Rational(3,10))},
              {0:alpha,5:s.sqrt(3)/5,10:-s.Rational(7,10)},
              {1:-s.Rational(7,10),6:s.sqrt(3)/5,11:alpha},
              {2:s.sqrt(s.Rational(3,10)),7:s.sqrt(s.Rational(7,10))}]
    V=s.Matrix(12,4,lambda r,a:supports[a].get(r,0))
    simplify=lambda A:A.applyfunc(lambda x:s.simplify(s.expand_complex(x)))
    assert s.simplify(alpha*s.conjugate(alpha))==s.Rational(39,100)
    assert simplify(V.H*V)==s.eye(4)
    splits=[]
    for a in range(4):
        splits.append(s.Matrix(2,11,lambda r,b:V[r+b,a]*s.sqrt(s.binomial(10,b)/s.binomial(11,r+b))))
    for a in range(4):
        for b in range(4):
            target=s.eye(2)/2 if a==b else s.zeros(2)
            assert simplify(splits[a]*splits[b].H)==target
    Jz=s.diag(*[s.Rational(11,2)-r for r in range(12)])
    Jplus=s.zeros(12)
    for r in range(1,12):Jplus[r-1,r]=s.sqrt(r*(12-r))
    for T in [Jz,Jplus,Jplus.T]:assert simplify(V.H*T*V)==s.zeros(4)
    # Verify that the simplified printed coefficient equals the supplied source.
    source=s.sqrt(39)*(12*s.sqrt(429)+5*s.I*s.sqrt(37609))/10010
    assert s.simplify(s.expand_complex(alpha-source))==0
    return {'passed':True,'n':11,'order':1,'dimension':4,
            'exactOrderedCrossReductionEntries':64,'exactSpinCompressionEntries':48,
            'orthonormalityEntries':16,'displayedCoefficientMatchesJuliaSource':True,
            'method':'Exact normalized-Dicke splitting and full cross reductions; independent compressed spin check.'}

if __name__=='__main__':
    result=verify()
    Path(__file__).with_name('lower-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
