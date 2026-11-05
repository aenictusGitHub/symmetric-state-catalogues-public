#!/usr/bin/env python3
"""Independent coefficientwise verifier: no irrep decomposition or numerical fit.

Uses a weighted binary-form basis e_r=sqrt(C(n,r))|D_n^r>.
For psi=sum z_r e_r, the output checks the whole polynomial identity
sum alpha_lmk f_lmk(psi) = (sum C(n,r)|z_r|^2)^3,
treating z and conjugate(z) as independent formal variables.
"""
import argparse
import itertools
import json
import math
import time
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import sympy as s
from sympy.physics.wigner import wigner_3j


def rational(x):
    x = s.cancel(x)
    assert x.is_Rational
    return Fraction(int(x.p), int(x.q))


def tensors(n):
    d = n + 1
    plus, minus = s.zeros(d), s.zeros(d)
    for r in range(n):
        plus[r+1,r] = r+1
        minus[r,r+1] = n-r
    out = {}
    for ell in range(d):
        matrix = plus**ell/s.factorial(ell)
        out[ell,ell] = matrix
        for q in range(ell, -ell, -1):
            matrix = (minus*matrix-matrix*minus)/(ell-q+1)
            assert all(x.is_Integer for x in matrix)
            out[ell,q-1] = matrix
    return out


@lru_cache(None)
def couplings(ell,m,k):
    entries=[]
    for q in range(-ell,ell+1):
        for r in range(-m,m+1):
            z=-q-r
            if abs(z)>k:
                continue
            value=wigner_3j(ell,m,k,q,r,z)
            if value:
                value/=s.sqrt(s.binomial(2*ell,ell-q)*s.binomial(2*m,m-r)*s.binomial(2*k,k-z))
                entries.append((q,r,z,value))
    seed=entries[0][-1]
    return [(q,r,z,rational(s.simplify(v/seed))) for q,r,z,v in entries]


def verify(certificate):
    start=time.time()
    n=certificate['n'];t=certificate['t'];d=n+1
    T=tensors(n)
    bilinear={key:[(row,col,int(s.binomial(n,row)*matrix[row,col]))
                   for row in range(d) for col in range(d) if matrix[row,col]]
              for key,matrix in T.items()}
    polynomial=defaultdict(Fraction)
    contributions=0
    for term in certificate['terms']:
        ell,m,k=term['ranks'];alpha=Fraction(term['coefficient'])
        assert any(1<=a<=t for a in (ell,m,k))
        assert (ell+m+k)%2==0
        for q,r,z,coupling in couplings(ell,m,k):
            coefficient=alpha*coupling
            for a,b,c in itertools.product(bilinear[ell,q],bilinear[m,r],bilinear[k,z]):
                bar=tuple(sorted((a[0],b[0],c[0])))
                ket=tuple(sorted((a[1],b[1],c[1])))
                polynomial[bar,ket]+=coefficient*a[2]*b[2]*c[2]
                contributions+=1
    for indices in itertools.combinations_with_replacement(range(d),3):
        multiplicity=math.factorial(3)
        for index in set(indices):
            multiplicity//=math.factorial(indices.count(index))
        polynomial[indices,indices]-=multiplicity*math.prod(math.comb(n,r) for r in indices)
    coherent_sum=Fraction(0)
    for term in certificate['terms']:
        ell,m,k=term['ranks']
        coherent=next(v for q,r,z,v in couplings(ell,m,k) if q==r==z==0)*math.comb(n,ell)*math.comb(n,m)*math.comb(n,k)
        assert coherent==Fraction(term['coherentValue'])
        normalized=Fraction(term['coefficient'])*coherent
        assert normalized==Fraction(term['coherentNormalizedCoefficient'])
        coherent_sum+=normalized
    assert coherent_sum==1
    nonzero={str(key):str(value) for key,value in polynomial.items() if value}
    result={'n':n,'t':t,'method':'Exact coefficientwise formal polynomial expansion in weighted binary-form coordinates; dense commutator construction of all tensor operators; no highest-weight block assumptions.',
            'terms':len(certificate['terms']),'monomialSlots':len(polynomial),
            'expandedContributions':contributions,'nonzeroRemainders':len(nonzero),
            'verifiedExactly':not nonzero,'runtimeSeconds':time.time()-start}
    if nonzero:
        result['firstRemainders']=dict(list(nonzero.items())[:10])
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('certificate',type=Path,nargs='?',default=Path(__file__).with_name('certificate.json'))
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('verification.json'))
    args=parser.parse_args()
    result=verify(json.loads(args.certificate.read_text()))
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    assert result['verifiedExactly'], result
