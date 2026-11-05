"""General AC-subspace bounds, with exact covariant-channel equality tests.

The proof is in proof.tex. No numerical zero tests are used.
"""
from functools import lru_cache
import sympy as s
from sympy.physics.wigner import wigner_6j

def packing_bound(n,t):
    if not 0<=t<=n:raise ValueError('Require 0 <= t <= N.')
    if t==0:return n+1
    r=t//2
    if t%2:return n//((r+1)*(r+2))
    d=(r+1)**2
    return 0 if n<d else (n+1)//d

@lru_cache(None)
def channel_multiplier(n,t,L):
    """HS-normalized rectangular channel; only exact vanishing is needed."""
    j=s.Rational(n,2);r=t//2
    if t%2:
        jp=j-s.Rational(1,2)
        return s.simplify(sum((-1)**(n+L+q)*(2*q+2)*
            wigner_6j(j,j,L,jp,jp,s.Rational(2*q+1,2)) for q in range(r+1)))
    return s.simplify(sum((-1)**(n+q+L)*(2*q+1)*
        wigner_6j(j,j,q,j,j,L) for q in range(r+1)))

def bound(n,t):
    upper=packing_bound(n,t)
    result={'packingUpper':upper,'upper':upper,'saturationTested':False,
            'saturationExcluded':False}
    if not t or not upper:return result
    r=t//2;odd=t%2==1;count=(r+1)*(r+2) if odd else (r+1)**2
    saturated=upper*count==(n if odd else n+1)
    if not saturated:return result
    roots=[L for L in range(t+1,n if odd else n+1)
           if channel_multiplier(n,t,L)==0]
    result.update(saturationTested=True,saturatingDimension=upper,channelRoots=roots)
    reason=None
    if not roots and (not odd or n>=7):
        reason='top-rank-projector-obstruction' if odd else 'no-nonscalar-channel-root'
    if odd and upper==1 and r>=1:
        value=channel_multiplier(n,t,t+1)
        result['nextRankMultiplier']=str(value)
        if value!=0:reason='pure-state-next-order-obstruction'
    if reason:
        result.update(upper=upper-1,saturationExcluded=True,exclusionReason=reason)
    return result
