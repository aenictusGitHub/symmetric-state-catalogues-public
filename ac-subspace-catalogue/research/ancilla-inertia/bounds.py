"""Exact second-order AC-subspace bound from a spin-3/2 ancillary system.

For every N >= 8 and t >= 2: 4k <= N-2.
The all-size proof is in proof.tex; no numerical solver is involved.
"""
from fractions import Fraction

def upper_bound(n, t=2):
    """Return this theorem's integer upper bound, or None outside its domain."""
    if not isinstance(n,int) or not isinstance(t,int):
        raise TypeError('N and t must be integers.')
    if n<0 or not 0<=t<=n:
        raise ValueError('Require N >= 0 and 0 <= t <= N.')
    return (n-2)//4 if n>=8 and t>=2 else None

def eigenvalues(n):
    """Exact eigenvalues and multiplicities of the canonical quadratic witness.

The four-sector decomposition applies for N >= 3. At N=7 two sectors
have zero eigenvalue, so the strict-inertia theorem does not apply.
"""
    if not isinstance(n,int) or n<3:
        raise ValueError('The four-sector formula requires integer N >= 3.')
    return [
        (Fraction((n+2)*(5*n-3),8),n-2),
        (Fraction(-n*(n-7),8),n),
        (Fraction(-(3*n*n+n+6),8),n+2),
        (Fraction(-n*(n-7),8),n+4),
    ]

def bound(n,t=2):
    upper=upper_bound(n,t)
    result={'applies':upper is not None,'upper':upper,'inequality':'4k <= N-2',
            'domain':'N >= 8 and t >= 2','proofId':'second-order-ancilla-inertia'}
    if upper is not None:
        result['inertia']={'positive':n-2,'negative':3*n+6,'zero':0}
    return result
