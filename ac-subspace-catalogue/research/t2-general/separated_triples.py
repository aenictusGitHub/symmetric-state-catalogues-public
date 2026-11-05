#!/usr/bin/env python3
"""Exact separated-triple construction of second-order anticoherent subspaces.

For N >= 6 this produces floor((N+3)/9) orthonormal vectors.  It is a
lower-bound construction, not a claim of universal optimality.
All checks use rational arithmetic; the amplitudes are square roots of the
listed rational weights.  Reproduce the JSON artifacts with:
    python3 research/t2-general/separated_triples.py
"""
from fractions import Fraction as Q
from pathlib import Path
import json


def construction(n):
    if not isinstance(n, int) or n < 6:
        raise ValueError('The construction requires an integer N >= 6.')
    k = (n + 3) // 9
    h = (n - 3 * (k - 1)) // 2
    j = Q(n, 2)
    variance = j * (j + 1) / 3
    rows = []
    for i in range(k):
        indices = [3 * i + ell * h for ell in range(3)]
        a, b, c = [j - r for r in indices]
        weights = [(variance + b*c)/(2*h*h),
                   -(variance + a*c)/(h*h),
                   (variance + a*b)/(2*h*h)]
        rows.append({'indices': indices, 'weights': weights})
    return {'n': n, 't': 2, 'k': k, 'h': h, 'rows': rows}


def verify(result):
    n, k, h = (result[key] for key in ('n', 'k', 'h'))
    j = Q(n, 2)
    variance = j * (j + 1) / 3
    supports = []
    for row in result['rows']:
        indices, weights = row['indices'], row['weights']
        assert all(0 <= r <= n for r in indices)
        assert all(w > 0 for w in weights)
        assert sum(weights) == 1
        assert sum(w*(j-r) for w, r in zip(weights, indices)) == 0
        assert sum(w*(j-r)**2 for w, r in zip(weights, indices)) == variance
        supports.extend(indices)
    supports.sort()
    assert len(supports) == len(set(supports)) == 3*k
    assert min(b-a for a, b in zip(supports, supports[1:])) >= 3
    assert h >= 3*k
    assert n-2*h-3*(k-1) in (0, 1)
    # The gap check annihilates all compressed operators shifting Dicke
    # excitation by one or two; the moment checks annihilate rank 1,2
    # diagonal tensors.  Consequently every cross reduction is delta_ab I3/3.
    return True


def output(result):
    return {**result, 'rows': [
        {**row, 'weights': [str(w) for w in row['weights']],
         'amplitudes': ['sqrt(' + str(w) + ')' for w in row['weights']]}
        for row in result['rows']]}


if __name__ == '__main__':
    directory = Path(__file__).resolve().parent
    checked = list(range(6, 1001))
    for n in checked:
        verify(construction(n))
    examples = [output(construction(n)) for n in (35, 42)]
    (directory / 'separated-triples-bases.json').write_text(
        json.dumps({'construction': 'separated triples',
                    'dimension': 'floor((N+3)/9)',
                    'domain': 'integer N >= 6',
                    'status': 'exact analytical lower bound',
                    'examples': examples}, indent=2) + '\n')
    (directory / 'separated-triples-verification.json').write_text(
        json.dumps({'arithmetic': 'exact rational', 'verified_n_min': 6,
                    'verified_n_max': 1000, 'verified_count': len(checked),
                    'passed': True,
                    'checks': ['strictly positive rational weights',
                               'unit normalization', 'zero first spin moment',
                               'isotropic second spin moment',
                               'disjoint supports with global gap at least 3'],
                    'qualification': 'Finite exact checks supplement the all-N proof.'},
                   indent=2) + '\n')
    print('Verified exact conditions for every N=6,...,1000.')
    for example in examples:
        print(json.dumps(example, indent=2))
