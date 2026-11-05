#!/usr/bin/env python3
"""Exactly verify the radical Dicke states shipped with this catalogue.

The separately certified algebraic N=11 double-point state is checked by the
verifiers in n11-g2/. This script does not use decimal amplitudes as evidence.
Requires Python 3 and SymPy; run from any working directory.
"""
import argparse
import json
import math
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent


def verify(state):
    n = int(state['n'])
    order = int(state['achievedOrder'])
    coefficients = [s.sympify(value) for value in state['exactCoefficients']]
    assert len(coefficients) == n + 1, state['id']
    assert s.simplify(sum(c * s.conjugate(c) for c in coefficients) - 1) == 0, state['id']
    reductions = []
    for t in range(1, order + 1):
        schmidt = s.zeros(t + 1, n - t + 1)
        for a in range(t + 1):
            for b in range(n - t + 1):
                factor = s.sqrt(s.Rational(math.comb(t, a) * math.comb(n - t, b), math.comb(n, a + b)))
                schmidt[a, b] = coefficients[a + b] * factor
        residual = (schmidt * schmidt.conjugate().T - s.eye(t + 1) / s.Integer(t + 1)).applyfunc(s.simplify)
        assert residual == s.zeros(t + 1), (state['id'], t, residual)
        reductions.append({'t': t, 'identityVerifiedExactly': True, 'purity': str(s.Rational(1, t + 1))})
    return {'id': state['id'], 'n': n, 'achievedOrder': order, 'normalizationVerifiedExactly': True, 'reductions': reductions}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=HERE / 'states-input.json')
    parser.add_argument('--output', type=Path, default=HERE / 'states-verification.json')
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    states = source['states'] if isinstance(source, dict) else source
    assert len({state['id'] for state in states}) == len(states)
    results = [verify(state) for state in states]
    output = {
        'method': 'Exact symbolic normalization and normalized-Dicke Schmidt reductions; every reduced matrix through the claimed achieved order equals I/(t+1).',
        'stateCount': len(results),
        'reductionCount': sum(len(state['reductions']) for state in results),
        'allVerifiedExactly': True,
        'states': results,
    }
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({key: value for key, value in output.items() if key != 'states'}, indent=2))


if __name__ == '__main__':
    main()
