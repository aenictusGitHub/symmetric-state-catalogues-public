#!/usr/bin/env python3
"""Calculate the exact radical four-band basis without numerical optimization.

Usage:
  python3 construct_exact.py --n 40 --exact-json four-band-n40-exact.json

Python API:
  from construct_exact import exact_basis
  B = exact_basis(40)  # SymPy sparse matrix; columns are Dicke-basis vectors.

The JSON uses ordered named definitions to avoid repeatedly expanding phase
products. Every amplitude can be calculated exactly by evaluating these
expressions in order. Only rationals, square roots, products and I occur.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sympy as S
from four_bands import symbolic_weights


def rational_squared_diagonal(squared_lengths):
    """Find a rational diagonal squared using only exact rational comparisons.

Each pair permits D in the interval centered at A+B with radius 2 sqrt(AB).
The universal polynomial theorem guarantees a nonempty interior intersection.
Bisection between the rational centers must therefore terminate.
"""
    A, B, C, E = squared_lengths
    centers = [A+B, C+E]
    radii_squared = [4*A*B, 4*C*E]
    def triangles(D):
        return [R-(D-c)**2 for c,R in zip(centers,radii_squared)]
    for D in centers:
        if D>0 and all(t>=0 for t in triangles(D)):
            return D
    low, high = min(centers), max(centers)
    while True:
        D = (low+high)/2
        tests = triangles(D)
        if all(t>=0 for t in tests):
            return D
        # Move toward the center of a failed triangle interval. There cannot
        # be conflicting failed intervals, since their interiors intersect.
        failed = next(i for i,t in enumerate(tests) if t<0)
        if D<centers[failed]:
            low=D
        else:
            high=D


def exact_construction(n):
    if not isinstance(n,int) or n<40:
        raise ValueError('The four-band theorem requires an integer N >= 40.')
    q,remainder=divmod(n,8)
    k,h,start=q-3,2*q,4+remainder//2
    j=S.Rational(n,2)
    v=j*(j+1)/(3*h*h)
    supports=[[start+2*i+ell*h for ell in range(4)] for i in range(k)]
    weights=[[S.factor(w) for w in symbolic_weights((j-start-2*i)/h,v)]
             for i in range(k)]
    definitions=[]
    phases=[[S.Integer(1)]*4]
    phase_names=[[S.Symbol(f'p_{i}_{ell}') for ell in range(4)] for i in range(k)]
    for ell in range(4):
        definitions.append((phase_names[0][ell],S.Integer(1)))
    diagonals=[]
    edge_data=[]
    for i in range(k-1):
        squared=[]
        for ell in range(4):
            rr=supports[i][ell]
            squared.append(weights[i][ell]*weights[i+1][ell]*
                           (rr+1)*(rr+2)*(n-rr)*(n-rr-1))
        D=rational_squared_diagonal(squared)
        diagonals.append(D)
        increments=[]
        tests=[]
        for pair_index in range(2):
            A,B=squared[2*pair_index:2*pair_index+2]
            T=4*A*B-(D-A-B)**2
            assert T>=0
            U,V=D+A-B,D+B-A
            assert U*U+T==4*D*A and V*V+T==4*D*B
            assert U+V==2*D
            sign=1 if pair_index==0 else -1
            increments.extend([sign*(U+S.I*S.sqrt(T))/(2*S.sqrt(D*A)),
                               sign*(V-S.I*S.sqrt(T))/(2*S.sqrt(D*B))])
            tests.append(T)
        next_phases=[]
        for ell in range(4):
            zname=S.Symbol(f'z_{i}_{ell}')
            definitions.append((zname,increments[ell]))
            definitions.append((phase_names[i+1][ell],phase_names[i][ell]*zname))
            next_phases.append(S.Mul(phases[i][ell],increments[ell],evaluate=False))
        phases.append(next_phases)
        edge_data.append({'squared_lengths':squared,'squared_diagonal':D,
                          'triangle_discriminants':tests,'increments':increments})
    entries={(supports[i][ell],i):S.Mul(S.sqrt(weights[i][ell]),phases[i][ell],evaluate=False)
             for i in range(k) for ell in range(4)}
    basis=S.ImmutableSparseMatrix(n+1,k,entries)
    # Exact all-column normalization and rank-1/rank-2 diagonal compressions.
    for indices,ww in zip(supports,weights):
        assert all(w>0 for w in ww)
        assert sum(ww)==1
        assert sum(w*(j-r) for w,r in zip(ww,indices))==0
        assert sum(w*(j-r)**2 for w,r in zip(ww,indices))==j*(j+1)/3
    # The supports have one parity, disjoint columns and no cross-band shifts
    # by 2. Each surviving within-band adjacent compression is cancelled by
    # the two exact triangle identities checked above.
    assert len({r for row in supports for r in row})==4*k
    assert h-2*(k-1)==8
    return {'n':n,'k':k,'supports':supports,'weights':weights,'phases':phases,
            'phase_names':phase_names,'definitions':definitions,'edges':edge_data,
            'basis':basis}


def exact_basis(n):
    """Return the exact (N+1)-by-k SymPy sparse Dicke coefficient matrix."""
    return exact_construction(n)['basis']


def exact_json(construction):
    c=construction
    return {'n':c['n'],'t':2,'k':c['k'],
            'status':'exact analytical construction; optimality not asserted',
            'basis_convention':'Column i is u_i; Dicke row r is D_N^(r).',
            'expression_convention':'SymPy syntax; I is the imaginary unit. Evaluate ordered definitions before amplitudes.',
            'definitions':[{'name':str(name),'expression':S.sstr(expr)}
                           for name,expr in c['definitions']],
            'basis':[{'column':i,'components':[
                {'dicke_index':c['supports'][i][ell],
                 'squared_magnitude':str(c['weights'][i][ell]),
                 'amplitude':S.sstr(S.sqrt(c['weights'][i][ell])*c['phase_names'][i][ell])}
                for ell in range(4)]} for i in range(c['k'])],
            'edges':[{'from_column':i,'to_column':i+1,
                      'squared_lengths':[str(a) for a in edge['squared_lengths']],
                      'squared_diagonal':str(edge['squared_diagonal']),
                      'triangle_discriminants':[str(t) for t in edge['triangle_discriminants']]}
                     for i,edge in enumerate(c['edges'])],
            'exact_checks':{'positive_rational_weights':True,'normalization':True,
                            'first_and_second_moments':True,
                            'phase_norms_and_adjacent_closures':True}}


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,required=True)
    parser.add_argument('--exact-json',type=Path)
    args=parser.parse_args()
    construction=exact_construction(args.n)
    payload=exact_json(construction)
    if args.exact_json:
        args.exact_json.write_text(json.dumps(payload,indent=2)+'\n')
        print(f'Wrote an exact {construction["k"]}-dimensional AC2 basis for N={args.n} to {args.exact_json}.')
    else:
        print(json.dumps(payload,indent=2))
