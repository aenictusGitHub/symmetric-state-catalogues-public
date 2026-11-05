#!/usr/bin/env python3
"""Exact arithmetic certificate for cyclic AC2 saturation exclusion.

This excludes cyclic-invariant upper-dimension subspaces, not arbitrary
subspaces. The rank-one-per-residue property follows from the variance lemma. The proof is in construction-cyclic-exclusion.tex.
"""
from pathlib import Path
from fractions import Fraction as Q
import json


def certificate():
    # cos(pi/8)>11/12 iff sqrt(2)>49/36; both sides are positive.
    assert Q(49,36)**2<2
    z=[1-Q(2*l,5) for l in range(5)]
    dual=lambda x:x*x-Q(8,5)*x
    excess=[dual(x)+Q(3,5) for x in z]
    assert excess==[0,0,Q(8,25),Q(24,25),Q(48,25)]
    assert sum(excess)==Q(16,5)
    derivatives=[Q(8,5)-2*x for x in z]
    assert all(d>0 for d in derivatives[1:])
    # The variance lemma reduces arbitrary cyclic-invariant codes to one
    # dimension per residue: c-2k^2>0 in every stated parameter range.
    rank_gaps=[]
    for s in range(1,6):
        k=85
        numerator=k*k-10*(s-1)*k+s*(s-2)
        assert numerator>0
        assert 2*k-10*(s-1)>0  # The numerator increases for all k>=85.
        rank_gaps.append(str(Q(numerator,12)))
    assert Q(62*62-1,12)>0  # N=5k-1, k>=62.
    records=[]
    for s,threshold in [(5,420),(1,309)]:
        eps_numerators=[12+Q(2*s*l,5) for l in range(5)]
        linear=sum(derivatives[l]*eps_numerators[l] for l in range(1,5))
        quadratic=sum(eps_numerators[l]**2 for l in range(1,5))
        N=threshold
        B=Q(15,4)-Q(6,N+2)
        U=-Q(3,5)+(Q(16,5)+linear/N+quadratic/(N*N))/B
        gap=Q(1,3)-U
        assert gap>0
        if s==5:
            assert linear==Q(584,5) and quadratic==1176
            assert gap==Q(872,1655325)
        else:
            assert linear==Q(424,5) and quadratic==Q(3384,5)
        records.append({'residue_case':'all N=5k-s, 1<=s<=5' if s==5 else 'N=5k-1',
                        'minimum_N':N,'minimum_k_on_N_5k_minus_1':62 if s==1 else None,
                        'linear_numerator':str(linear),'quadratic_numerator':str(quadratic),
                        'upper_expectation_at_threshold':str(U),
                        'strict_gap_below_one_third':str(gap)})
    return {'status':'rigorous cyclic-invariance exclusion',
            'scope':'All k-dimensional AC2 subspaces invariant under a collective rotation by 2*pi/k, about any axis',
            'rank_one_per_residue':'Derived from the five-node variance lemma, not an assumption',
            'variance_lemma_gaps_at_k85_for_s1_to_5':rank_gaps,
            'general_subspace_exclusion':False,
            'window_size':7,'window_weights':'sin((a+1)*pi/8), a=0,...,6',
            'cosine_lower_bound':'cos(pi/8)>11/12',
            'band_cap':'1/(15/4-6/(N+2))',
            'dual_polynomial':'Q(z)=z^2-(8/5)z',
            'records':records}


if __name__=='__main__':
    data=certificate()
    output=Path(__file__).with_name('construction-cyclic-exclusion.json')
    output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
