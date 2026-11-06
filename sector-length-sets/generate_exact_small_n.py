#!/usr/bin/env python3
"""Analytical Dicke expressions matching the exported N=3,4,5 witnesses."""
from pathlib import Path
from fractions import Fraction
from math import sqrt, sin, cos, pi
import json
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
states = json.loads((DATA / 'small-n.json').read_text())


def D(n, k):
    return rf'\lvert D_{n}^{{({k})}}\rangle'


def rational(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return rf'\frac{{{q.numerator}}}{{{q.denominator}}}'


def pi_multiple(q):
    q = Fraction(q)
    if not q: return '0'
    if q == 1: return r'\pi'
    if q.numerator == 1: return rf'\frac{{\pi}}{{{q.denominator}}}'
    return rf'\frac{{{q.numerator}\pi}}{{{q.denominator}}}'


def interpolation(a, b, q):
    q = Fraction(q)
    if not q: return a
    if q == 1: return b
    return rf'{rational(1-q)}\left({a}\right)+{rational(q)}\left({b}\right)'


named = {}
for n in [3, 4, 5]:
    named[f'n{n}-product'] = {'latex': rf'\lvert0\rangle^{{\otimes{n}}}={D(n,0)}'}
    named[f'n{n}-ghz'] = {'latex': rf'\frac{{{D(n,0)}+{D(n,n)}}}{{\sqrt2}}'}
    named[f'n{n}-w'] = {'latex': D(n,1)}
for n in [4,5]:
    named[f'n{n}-dicke2'] = {'latex': D(n,2)}
named['n4-tetra'] = {'latex': rf'\frac{{{D(4,0)}+\sqrt2\,{D(4,3)}}}{{\sqrt3}}'}
named['n5-psia'] = {'latex': rf'\frac{{\sqrt5\,{D(5,1)}+{D(5,5)}}}{{\sqrt6}}'}
named['n5-psib'] = {'latex': rf'\frac{{\sqrt5\,{D(5,0)}+\sqrt2\,{D(5,2)}-3{D(5,4)}}}{{4}}',
                    'note':'The stored numerical representative differs by a global minus sign.'}
named['n5-psic'] = {'latex': rf'\frac{{{D(5,1)}+\lambda_c{D(5,5)}}}{{\sqrt{{1+\lambda_c^2}}}}',
                    'parameterLatex': r'\lambda_c=2-\frac3{\sqrt5}'}

rho_latex = r'\rho(p)=\frac{2p(5p^2+6\sqrt5\,p-11)}{3\sqrt5\,p^2+2p+3\sqrt5}'
ghz = lambda n: rf'\cos\frac\phi2\,{D(n,0)}+\sin\frac\phi2\,{D(n,n)}'
psi4 = rf'\cos\frac\theta2\!\left(\cos\frac\phi2\,{D(4,0)}+\sin\frac\phi2\,{D(4,4)}\right)+i\sin\frac\theta2\,{D(4,2)}'
left5 = rf'\frac{{{D(5,0)}+\sqrt{{3t^2-5}}\,{D(5,2)}+t{D(5,4)}}}{{2\sqrt{{t^2-1}}}}'
families = {
    'n3-entire-interval': {'latex':ghz(3),'parameterLatex':r'0\le\phi\le\pi/2',
                           'sectorLatex':r'S_1=3\cos^2\phi,\quad S_2=3,\quad S_3=4-S_1'},
    'n4-lower': {'latex':psi4,
                 'parameterLatex':r'u=x/4,\quad\theta=\arccos u,\quad\phi=\arcsin\!\frac{1-u}{1+u}',
                 'sectorLatex':r'S_1=x,\quad S_2=x+2,\quad0\le x\le4'},
    'n4-left': {'latex':psi4,
                'parameterLatex':r'\theta=\arccos\sqrt{(y-2)/4},\quad\phi=\pi/2',
                'sectorLatex':r'S_1=0,\quad S_2=y,\quad2\le y\le6'},
    'n4-top': {'latex':ghz(4),'parameterLatex':r'0\le\phi\le\pi/2',
               'sectorLatex':r'S_1=4\cos^2\phi,\quad S_2=6'},
    'n4-exact-region': {'latex':psi4,
                        'parameterLatex':r'Y=3-y/2,\quad u=\frac{x+Y\sqrt{4-x-2Y}}{x+2Y},\quad v=1-\frac{x(1-u)}{Y(1+u)},\quad\theta=\arccos u,\quad\phi=\arcsin v',
                        'note':'For Y=0 use u=1 and v=√(1−x/4). This inverse covers the entire exact triangle.'},
    'n5-small-arc': {'latex':rf'\frac{{-p{D(5,0)}+\sqrt{{\rho(p)}}\,{D(5,2)}+{D(5,4)}}}{{\sqrt{{1+p^2+\rho(p)}}}}',
                     'parameterLatex':rho_latex,
                     'domainLatex':r'2-3/\sqrt5\le p\le\sqrt5/3'},
    'n5-root-arc': {'latex':rf'\frac{{{D(5,1)}+\lambda{D(5,5)}}}{{\sqrt{{1+\lambda^2}}}}',
                    'parameterLatex':r'0\le\lambda\le2-3/\sqrt5',
                    'sectorLatex':r'S_1=\frac{(3-5\lambda^2)^2}{5(1+\lambda^2)^2},\quad S_2=\frac{2(25\lambda^4+10\lambda^2+9)}{5(1+\lambda^2)^2}'},
    'n5-linear-arc': {'latex':rf'\cos\theta\,{D(5,0)}+\sin\theta\,{D(5,1)}',
                      'parameterLatex':r'0\le\theta\le\pi/2',
                      'sectorLatex':r'S_1=5-\frac{16}{5}\sin^4\theta,\quad S_2=2S_1'},
    'n5-top': {'latex':ghz(5),'parameterLatex':r'0\le\phi\le\pi/2',
               'sectorLatex':r'S_1=5\cos^2\phi,\quad S_2=10'},
    'n5-left-negative': {'latex':left5,
                         'parameterLatex':r'-\sqrt{(345+60\sqrt{15})/289}\le t\le-3/\sqrt5',
                         'sectorLatex':r'S_1=0,\quad18/5\le S_2\le15/4'},
    'n5-left-positive': {'latex':left5,
                         'parameterLatex':r'\sqrt{5/3}\le t\le\sqrt5',
                         'sectorLatex':r'S_1=0,\quad15/4\le S_2\le10'},
}

# Exact grid values, rather than equating a rounded floating-point value to a
# symbolic parameter. The float values remain available in the witness record.
rt5 = sqrt(5)
pc = 2-3/rt5
p0 = rt5/3
checks = []
for nstr, d in states.items():
    n = int(nstr)
    for curve in d['curves']:
        cid = curve['id']
        family = families[cid]
        size = len(curve['points'])-1
        for j, witness in enumerate(curve['points']):
            q = Fraction(j,size)
            entry = {'latex':family['latex'],'familyId':cid,'analyticalStatus':'exact state family and exact generator-grid parameter'}
            if cid.endswith('-top') or cid=='n3-entire-interval':
                entry['parameterLatex'] = rf'\phi={pi_multiple(q/2)}'
                value = float(q)*pi/2
                param = 'phi'
            elif cid == 'n4-lower':
                entry['parameterLatex'] = rf'x={rational(4*q)},\quad u=x/4,\quad\theta=\arccos u,\quad\phi=\arcsin\frac{{1-u}}{{1+u}}'
                value = float(4*q)
                param = 'S1'
            elif cid == 'n4-left':
                entry['parameterLatex'] = rf'y={rational(2+4*q)},\quad\theta=\arccos\sqrt{{(y-2)/4}},\quad\phi=\pi/2'
                value = float(2+4*q)
                param = 'S2'
            elif cid == 'n5-small-arc':
                exactp = interpolation(r'\sqrt5/3',r'2-3/\sqrt5',q)
                entry['parameterLatex'] = rf'p={exactp},\qquad {rho_latex}'
                value = (1-float(q))*p0+float(q)*pc
                param = 'p'
            elif cid == 'n5-root-arc':
                exactlam = '0' if q==1 else (r'2-3/\sqrt5' if q==0 else rf'{rational(1-q)}\left(2-3/\sqrt5\right)')
                entry['parameterLatex'] = rf'\lambda={exactlam}'
                value = (1-float(q))*pc
                param = 'lambda'
            elif cid == 'n5-linear-arc':
                entry['parameterLatex'] = rf'\theta={pi_multiple((1-q)/2)}'
                value = float(1-q)*pi/2
                param = 'theta'
            elif cid == 'n5-left-negative':
                exactt = interpolation(r'-3/\sqrt5',r'-\sqrt{(345+60\sqrt{15})/289}',q)
                entry['parameterLatex'] = rf't={exactt}'
                value = -(1-float(q))*3/rt5-float(q)*sqrt((345+60*sqrt(15))/289)
                param = 't'
            elif cid == 'n5-left-positive':
                exactt = interpolation(r'\sqrt{5/3}',r'\sqrt5',q)
                entry['parameterLatex'] = rf't={exactt}'
                value = (1-float(q))*sqrt(5/3)+float(q)*rt5
                param = 't'
            else:
                raise ValueError(cid)
            if 'sectorLatex' in family:
                entry['sectorLatex'] = family['sectorLatex']
            named[witness['id']] = entry
            checks.append(abs(value-witness['parameter'][param]))

assert max(checks)<2e-14
payload=json.dumps(named,ensure_ascii=False,separators=(',',':'))
family_payload=json.dumps(families,ensure_ascii=False,separators=(',',':'))
(DATA/'exact-small-n.js').write_text('/* Analytical expressions for the actual stored Dicke representatives, up to global phase only. */\n'
    +'globalThis.EXACT_SMALL_EXPRESSIONS = '+payload+';\n'
    +'globalThis.EXACT_BOUNDARY_FAMILIES = '+family_payload+';\n')
(DATA/'exact-small-n.json').write_text(json.dumps({'states':named,'families':families},ensure_ascii=False,indent=2)+'\n')

# Check the named formula vectors against the actual stored representatives,
# allowing only global phase (never a change of collective frame).
expected={}
def vec(n,terms):
    c=np.zeros(n+1,dtype=complex)
    for k,v in terms.items():c[k]=v
    return c/np.linalg.norm(c)
for n in [3,4,5]:
    expected[f'n{n}-product']=vec(n,{0:1})
    expected[f'n{n}-ghz']=vec(n,{0:1,n:1})
    expected[f'n{n}-w']=vec(n,{1:1})
for n in [4,5]:expected[f'n{n}-dicke2']=vec(n,{2:1})
expected['n4-tetra']=vec(4,{0:1,3:sqrt(2)})
expected['n5-psia']=vec(5,{1:sqrt(5),5:1})
expected['n5-psib']=vec(5,{0:sqrt(5),2:sqrt(2),4:-3})
expected['n5-psic']=vec(5,{1:1,5:pc})
max_state_error=0
for d in states.values():
    for r in d['knownStates']:
        actual=np.array([complex(*z)for z in r['coefficients']]);actual/=np.linalg.norm(actual)
        v=expected[r['id']]
        phase=np.vdot(v,actual)
        error=np.max(np.abs(actual-v*phase/abs(phase)))
        max_state_error=max(max_state_error,float(error))
assert max_state_error<2e-12
audit={'status':'PASS','namedStates':len(expected),'familyExpressions':len(families),
       'totalPointExpressions':len(named),'exactGridParametersChecked':len(checks),
       'maxGridParameterError':max(checks),'maxNamedStateErrorUpToGlobalPhase':max_state_error,
       'note':'All exact sample parameters follow the analytical sampling grid; no decimal is presented as an exact parameter.'}
(DATA/'exact-small-n-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
