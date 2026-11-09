"""Build the two exceptional quintic covariants from the cubic purity gradient."""
from collections import defaultdict
from math import comb
from fractions import Fraction

def gradient_highest(J, beta):
    """Return the highest-weight polynomial with exact coefficients for a given beta."""
    if J not in (31, 33):
        raise ValueError('Expected doubled spin 31 or 33.')
    gradient = []
    for k in range(2):
        p = defaultdict(lambda: 0)
        for i in range(5):
            for h in range(5):
                for j in range(8):
                    for l in range(8):
                        if h+j == k:
                            a = tuple(sorted((i+j,h+l)))
                            p[a,(i+l,)] += 2*comb(4,i)*comb(4,h)*comb(7,j)*comb(7,l)
        for m in range(12):
            a = tuple(sorted((k,m)))
            p[a,(m,)] -= 2*beta*comb(11,k)*comb(11,m)
        gradient.append({(a,b):v*Fraction(1,comb(11,k)) for (a,b),v in p.items() if v})
    terms = [(0,(0,0),1)] if J==33 else [(1,(0,0),1),(0,(0,1),-1)]
    out = defaultdict(lambda: 0)
    for k,prefix,sign in terms:
        for (a,b),v in gradient[k].items():
            out[tuple(sorted(prefix+a)),b] += sign*v
    return {m:v for m,v in out.items() if v}
