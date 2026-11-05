"""Independently check exact S4 attainment by three source constructions.

The global lower bound S4 >= binomial(N,4)/5 is the cited multipole
identity. This script checks normalization and equality, using reduced
purities instead of the source repository's moment verifier.
"""
from math import comb
import sympy as s

def representative(n):
    c = [s.S.Zero] * (n + 1)
    if n == 10:
        u = (-179 + s.I * s.sqrt(42479)) / (18 * s.sqrt(230))
        v = (12 * s.sqrt(115) + 30 * s.sqrt(2) * u) / (40 * s.sqrt(7))
        c[0], c[3] = s.sqrt(s.Rational(23,126)), s.sqrt(s.Rational(3,7))
        c[6], c[9] = u / s.sqrt(6), s.expand(u * v * s.sqrt(2) / 3)
    elif n == 11:
        c[0], c[5], c[10] = map(s.sqrt, [s.Rational(7,150), s.Rational(11,25), s.Rational(77,150)])
    elif n == 12:
        c[0], c[5], c[10] = map(s.sqrt, [s.Rational(32,175), s.Rational(66,175), s.Rational(11,25)])
    else:
        raise ValueError(n)
    return c

def exact_sector(c, k):
    n = len(c)-1
    norm = s.simplify(sum(z*s.conjugate(z) for z in c))
    assert norm == 1
    purities = []
    for r in range(k+1):
        a = s.Matrix(r+1,n-r+1,lambda i,j: c[i+j] * s.sqrt(s.Rational(comb(r,i)*comb(n-r,j),comb(n,i+j))))
        rho = a*a.conjugate().T
        purities.append(s.simplify(s.trace(rho*rho)))
    value = s.simplify(comb(n,k) * sum((-1)**(k-r)*comb(k,r)*2**r*purities[r] for r in range(k+1)))
    return norm,value

if __name__ == '__main__':
    for n in (10,11,12):
        norm,value=exact_sector(representative(n),4)
        assert value == s.Rational(comb(n,4),5)
        print(f'PASS N={n}: exact norm={norm}, S4={value}=C({n},4)/5')
