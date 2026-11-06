#!/usr/bin/env python3
"""Exact metadata for the filtered-octa C3 rim and rank-four equality line.

This helper never edits the atlas, its existing data, or the manuscript.
All input parameters are rational by construction.  Sector values follow
the exact image identities proved with the original state constructors.
An independent reduced-density-matrix calculation audits each witness.
"""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
OUTPUT = ROOT / 'data/n6-boundary-sample-exact.json'
AUDIT = ROOT / 'audit/n6-boundary-sample-exact-audit.json'


def full_tuple(x, y, z):
    return [x, y, z, 45+10*x-2*y-3*z, 3*z-9*x, 18-2*x+y-z]


def latex(v):
    return sp.latex(sp.factor(v))


def independent_sectors(c):
    """Get all sectors by binomial inversion of reduced-state purities."""
    c = np.asarray(c, dtype=complex)
    out = []
    for k in range(1, 7):
        a = np.zeros((k+1, 7-k), dtype=complex)
        for i in range(k+1):
            for j in range(7-k):
                a[i, j] = c[i+j] * math.sqrt(
                    math.comb(k, i)*math.comb(6-k, j)/math.comb(6, i+j))
        rho = a @ a.conj().T
        purity = float(np.trace(rho @ rho).real)
        normalized = 2**k*purity-1
        for l, sector in enumerate(out, start=1):
            normalized -= math.comb(k, l)/math.comb(6, l)*sector
        out.append(math.comb(6, k)*normalized)
    return out


def record(sid, values, expression, numerical_coefficients, provenance):
    expected = [float(v) for v in values]
    actual = independent_sectors(numerical_coefficients)
    err = max(abs(a-b) for a, b in zip(expected, actual))
    norm_error = abs(sum(abs(z)**2 for z in numerical_coefficients)-1)
    assert err < 2e-10, (sid, err)
    assert norm_error < 2e-13, (sid, norm_error)
    exact[sid] = {
        'analyticalExpression': expression,
        'exactSectorValues': {
            'values': [latex(v) for v in values],
            'status': 'exact',
            'provenance': provenance,
        },
    }
    audit[sid] = {'normalizationError': norm_error,
                  'maximumSectorError': err,
                  'exactSectorValues': [str(sp.factor(v)) for v in values],
                  'numericOracle': actual}


exact, audit = {}, {}
fo_formula = (
    r'\begin{aligned}|F_t\rangle={}&'
    r'\sqrt{\frac{4+5t+3\sqrt{t(t+8)}}{18}}\,D_0\\'
    r'&+\frac{\sqrt{5(1-t)}}{3}\,D_3\\'
    r'&-\sqrt{\frac{4+5t-3\sqrt{t(t+8)}}{18}}\,D_6'
    r'\end{aligned}')
for j in range(1, 9):
    t = sp.Rational(j, 9)
    x = sp.Rational(2, 3)*t*(t+8)
    y = 5+10*t*t
    z = sp.Rational(20, 3)*t*(4-t)
    coeff = [0.]*7
    coeff[0] = float(sp.sqrt((4+5*t+3*sp.sqrt(t*(t+8)))/18))
    coeff[3] = float(sp.sqrt(5*(1-t))/3)
    coeff[6] = -float(sp.sqrt((4+5*t-3*sp.sqrt(t*(t+8)))/18))
    sid = f'n6-fo-c3-{j}'
    record(sid, full_tuple(x, y, z), {
        'latex': fo_formula,
        'parameterLatex': r't='+latex(t)+r',\quad D_k=|D_6^{(k)}\rangle',
        'note': 'Exact normalized C3 state on the filtered-octahedral boundary.',
        'proofReference': 'Julia/N6_filtered_octa_attainment.jl, n6fo_c3_state and n6fo_c3_point',
    }, coeff, 'Exact filtered-octahedral C3 image; rational t=j/9.')

family_source = WORKSPACE / 'output/N6_plane_equality_segment_20261021/real/fourier_equality_family.py'
spec = importlib.util.spec_from_file_location('rank4_exact_original', family_source)
family = importlib.util.module_from_spec(spec)
spec.loader.exec_module(family)

rank4_formula = (
    r'|R_t\rangle=\frac{1}{\sqrt{\mathcal N}}'
    r'\sum_{k=0}^{6}\frac{p_k}{\sqrt{\binom{6}{k}}}\,D_k')
rank4_defs = r'''\begin{aligned}
g&=13-7t^2,\\
\Delta&=427-639t+287t^2-35t^3,\\
a&=\sqrt{\frac{35(7-5t^2)}{16g}},\\
R&=\frac{25(23t^2-84t+91)}{108g},\\
U&=-\frac{3(7-3t)g}{\Delta},\quad
V=-\frac{24a(t-1)g}{5\Delta},\\
P&=\frac{5a(1985-3484t+1594t^2+84t^3-119t^4)}{54g^2},\\
d&=\sqrt{\frac{R+\sqrt{R^2-4P^2}}{2}},\quad e=\frac{P}{d},\\
b&=Ud+Ve,\quad c=-Vd-Ue,\\
\mathcal N&=\frac{64(5-2t^2)}{27g},\\
p_0&=\frac a3-\frac c4-\frac e{20},\\
p_1&=\frac b2-1+\frac{3d}{10},\\
p_2&=a-\frac c4+\frac{3e}4,\\
p_3&=b-d,\\
p_4&=a+\frac c4-\frac{3e}4,\\
p_5&=\frac b2+1+\frac{3d}{10},\\
p_6&=\frac a3+\frac c4+\frac e{20},\\
D_k&=|D_6^{(k)}\rangle.
\end{aligned}'''
params = [sp.Rational(41, 40), sp.Rational(21, 20), sp.Rational(43, 40),
          sp.Rational(11, 10), sp.Rational(9, 8), sp.Rational(23, 20),
          sp.Rational(47, 40)]
for j, t in enumerate(params, start=1):
    r = 13*(t*t-1)/(4*(5-2*t*t))
    q = r*r
    x = sp.Rational(3, 2)-sp.Rational(474, 169)*q
    y = sp.Rational(99, 10)+sp.Rational(3234, 845)*q
    z = 11-sp.Rational(2300, 169)*q
    coeff = family.state(float(t))
    sid = f'n6-rank4-segment-{j}'
    record(sid, full_tuple(x, y, z), {
        'latex': rank4_formula,
        'parameterLatex': r't='+latex(t)+r'\qquad '+rank4_defs,
        'note': 'Explicit normalized real-radical state on the attained equality segment.',
        'proofReference': 'Julia/N6_rankfour_equality_segment.jl, n6rank4_equality_state; exact Fourier-family multipole identities',
    }, coeff, 'Exact rank-four segment image with q=[13(t^2-1)/(4(5-2t^2))]^2.')
    exact[sid]['exactSectorValues']['definitionsLatex'] = (
        r't='+latex(t)+r',\quad q=\left['
        r'\frac{13(t^2-1)}{4(5-2t^2)}\right]^2='+latex(q))

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
AUDIT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(exact, indent=2, ensure_ascii=False)+'\n')
result = {'status': 'PASS', 'states': len(exact),
          'maximumSectorError': max(r['maximumSectorError'] for r in audit.values()),
          'entries': audit}
AUDIT.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'entries'}, indent=2))
