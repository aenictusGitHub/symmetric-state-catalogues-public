#!/usr/bin/env python3
"""Exact sector readouts derived from the recorded analytical state families.

No rational reconstruction of rounded sector data is used. Saved floating-point
sector values serve only as an independent final numerical check.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
R = sp.Rational
small = json.loads((DATA / 'small-n.json').read_text())
n6 = json.loads((DATA / 'n6.json').read_text())
small_exact = json.loads((DATA / 'exact-small-n.json').read_text())['states']
all_records = {p['id']: p for d in small.values() for p in d['knownStates']}
all_records.update({p['id']: p for p in n6['states']})
for d in small.values():
    for curve in d['curves']:
        all_records.update({p['id']: p for p in curve['points']})

entries = {}
audits = {}


def tex(expr):
    """Keep a rational/radical value in a single readable exact fraction."""
    expr = sp.together(sp.radsimp(sp.simplify(expr)))
    numerator, denominator = sp.fraction(expr)
    return (sp.latex(numerator) if denominator == 1 else
            r'\frac{' + sp.latex(numerator) + '}{' + sp.latex(denominator) + '}')


def add(sid, expressions, *, definitions='', provenance='', values=None,
        substitutions=None):
    record = all_records[sid]
    evaluated = [float(sp.N(e.subs(substitutions or {}), 35)) for e in expressions]
    errors = [abs(v - w) for v, w in zip(evaluated, record['sectors'])]
    assert len(expressions) == len(record['coefficients']) - 1
    assert max(errors) < 3e-9, (sid, errors, evaluated, record['sectors'])
    entry = {'values': values or [tex(e) for e in expressions],
             'status': 'exact', 'provenance': provenance}
    if definitions:
        entry['definitionsLatex'] = definitions
    entries[sid] = entry
    audits[sid] = {'maximumAbsoluteSectorError': max(errors),
                   'exactExpressionValues': evaluated,
                   'storedNumericalValues': record['sectors']}


def tuple5(x, y):
    return [x, y, 10 + 2*x - y, 15-y, 6-3*x+y]


def tuple6(x, y, z):
    return [x, y, z, 45+10*x-2*y-3*z, 3*z-9*x, 18-2*x+y-z]


# The previously recorded exact small-size tuples are source data, not inferred
# from their displayed decimal values.
for d in small.values():
    for record in d['knownStates']:
        if record.get('exactSectors'):
            exact = []
            for value in record['exactSectors']:
                q = Fraction(value)
                exact.append(R(q.numerator, q.denominator))
            add(record['id'], exact,
                provenance='Recorded analytical representative and exact sector tuple; verified by reduced purities.')

add('n5-dicke2', tuple5(R(1,5), R(38,5)),
    provenance='Exact Dicke-state reduced purities.')
lam_c = 2 - 3/sp.sqrt(5)
x_c = (3-5*lam_c**2)**2/(5*(1+lam_c**2)**2)
y_c = 2*(25*lam_c**4+10*lam_c**2+9)/(5*(1+lam_c**2)**2)
add('n5-psic', tuple5(x_c, y_c),
    provenance='Exact square-root boundary family at lambda_c=2-3/sqrt(5).')

# All seventeen manuscript N=6 representatives have closed exact values.
named6 = {
    'n6-product': (6,15,20),
    'n6-ghz': (0,15,0),
    'n6-w': (R(8,3),5,R(40,3)),
    'n6-octa': (0,5,0),
    'n6-dicke-2': (R(2,3),R(43,5),R(28,3)),
    'n6-c4': (0,R(27,5),8),
    'n6-omega4': (R(1,6),5,R(25,3)),
    'n6-psi6': (R(24,49),R(405,49),R(440,49)),
    'n6-chi6': (R(24,49),R(405,49),R(440,49)),
    'n6-chi7-10': (R(3,2),R(99,10),11),
    'n6-chi1-2': (R(1,6),R(15,2),R(25,3)),
    'n6-ac2': (0,5,R(640,81)),
    'n6-psi0': (0,R(51,10),8),
    'n6-xi': (0,13-R(8,5)*sp.sqrt(15),8),
    'n6-omega': (R(63,121),R(618,55),R(756,121)),
    'n6-eta6': (R(2,27),R(65,9),R(220,27)),
}
for sid, (x,y,z) in named6.items():
    add(sid, [sp.sympify(v) for v in tuple6(x,y,z)],
        provenance='Analytical normalized Dicke representative, spin moments/spin flip, and pure-state complementary-sector identities.')

# E_*: f_*^2 and q_*=|Theta|^2 are exact manuscript parameters.
# z=8+2x-8q and the sharp supporting-plane equality determine y. This
# simplifies the nested cofactor state into values in Q(sqrt(21)).
ell = sp.sqrt(21)
x = (-187677+44793*ell)/32250
q = R(14027,5375)-R(10647,21500)*ell
z = 8+2*x-8*q
y = (1011-63*z+217*x)/65
add('n6-estar', tuple6(x,y,z),
    provenance='Exact exceptional-orbit parameters f_* and |Theta|^2, plus 65S2+63S3-217S1=1011; all values simplify in Q(sqrt(21)).')


# Reproduce only the explicitly displayed forty-one N=5 boundary samples.
def near(a,b):
    return math.hypot(a['sectors'][0]-b['sectors'][0],
                      a['sectors'][1]-b['sectors'][1]) < .02


chosen = []
n5 = small['5']
for curve in n5['curves']:
    count = len(curve['points'])-1
    for j in range(1,10):
        idx = int(j*count/10+.5)
        record = curve['points'][idx]
        if any(near(record,p) for p in chosen+n5['knownStates']):
            continue
        chosen.append(record)
        cid, sid = curve['id'], record['id']
        grid = R(idx,count)
        metadata = small_exact[sid]
        parameter = metadata['parameterLatex']
        if cid == 'n5-top':
            phi = sp.Symbol('phi', real=True)
            x = 5*sp.cos(phi)**2
            expressions = tuple5(x, sp.Integer(10))
            values = [r'5\cos^2\phi', '10', r'10\cos^2\phi', '5', r'16-15\cos^2\phi']
            add(sid, expressions, values=values, definitions=parameter,
                substitutions={phi:grid*sp.pi/2},
                provenance='Generalized GHZ family; exact rational multiple of pi from the generator grid.')
        elif cid == 'n5-linear-arc':
            theta = sp.Symbol('theta', real=True)
            x = 5-R(16,5)*sp.sin(theta)**4
            values = [r'5-\frac{16}{5}\sin^4\theta', r'10-\frac{32}{5}\sin^4\theta',
                      '10', r'5+\frac{32}{5}\sin^4\theta', r'1+\frac{16}{5}\sin^4\theta']
            add(sid, tuple5(x,2*x), values=values, definitions=parameter,
                substitutions={theta:(1-grid)*sp.pi/2},
                provenance='Exact product–W boundary family and generator-grid angle.')
        elif cid == 'n5-root-arc':
            lam = sp.Symbol('lambda', real=True)
            x = (3-5*lam**2)**2/(5*(1+lam**2)**2)
            y = 2*(25*lam**4+10*lam**2+9)/(5*(1+lam**2)**2)
            add(sid, tuple5(x,y), definitions=parameter,
                substitutions={lam:(1-grid)*lam_c},
                provenance='Exact two-level C4 boundary family and generator-grid lambda.')
        elif cid == 'n5-small-arc':
            p = sp.Symbol('p', real=True)
            rho = 2*p*(5*p**2+6*sp.sqrt(5)*p-11)/(3*sp.sqrt(5)*p**2+2*p+3*sp.sqrt(5))
            # These moment formulas are evaluated before substituting rho.
            rhos = sp.Symbol('rho', real=True)
            den = 1+p**2+rhos
            x = (5*p**2+rhos-3)**2/(5*den**2)
            y = (2*(5*p**2-rhos+1)**2+(4+6*rhos)**2+8*rhos*(3-sp.sqrt(5)*p)**2)/(5*den**2)
            exactp = (1-grid)*sp.sqrt(5)/3+grid*lam_c
            definitions = parameter + r',\quad D=1+p^2+\rho(p),\quad ' + \
                r'\begin{aligned}x_p&=\frac{(5p^2+\rho(p)-3)^2}{5D^2},\\' + \
                r'y_p&=\frac{2(5p^2-\rho(p)+1)^2+(4+6\rho(p))^2+8\rho(p)(3-\sqrt5p)^2}{5D^2}.\end{aligned}'
            add(sid, tuple5(x,y), values=['x_p','y_p','10+2x_p-y_p','15-y_p','6-3x_p+y_p'],
                definitions=definitions, substitutions={p:exactp,rhos:rho.subs(p,exactp)},
                provenance='Exact small C2 arc; dipole and quadrupole formulas from its three normalized Dicke amplitudes.')
        elif cid.startswith('n5-left-'):
            t = sp.Symbol('t', real=True)
            y = (177*t**4+36*sp.sqrt(5)*t**3-410*t**2-60*sp.sqrt(5)*t+225)/(20*(t**2-1)**2)
            exactt = ((1-grid)*(-3/sp.sqrt(5))+grid*(-sp.sqrt((345+60*sp.sqrt(15))/289))
                      if cid.endswith('negative') else
                      (1-grid)*sp.sqrt(R(5,3))+grid*sp.sqrt(5))
            definitions = parameter + r',\qquad ' + \
                r'Y(t)=\frac{177t^4+36\sqrt5t^3-410t^2-60\sqrt5t+225}{20(t^2-1)^2}'
            add(sid, tuple5(sp.Integer(0),y), values=['0','Y(t)','10-Y(t)','15-Y(t)','6+Y(t)'],
                definitions=definitions, substitutions={t:exactt},
                provenance='Exact zero-dipole family; quadrupole formula and exact generator-grid t.')
        else:
            raise ValueError(cid)

assert len(chosen)==41
assert entries['n6-c4']['values']==['0',r'\frac{27}{5}','8',r'\frac{51}{5}','24',r'\frac{77}{5}']

payload = json.dumps(entries, ensure_ascii=False, separators=(',',':'))
(DATA/'exact-sector-values.js').write_text(
    '/* Exact sector values from analytical representatives; no decimal rational fitting. */\n'
    + 'globalThis.EXACT_SECTOR_VALUES='+payload+';\n')
(DATA/'exact-sector-values.json').write_text(json.dumps(entries,ensure_ascii=False,indent=2)+'\n')
report={'status':'PASS','entries':len(entries),'N6ManuscriptRepresentatives':17,
        'N5DisplayedBoundarySamples':len(chosen),
        'maximumAbsoluteSectorError':max(v['maximumAbsoluteSectorError'] for v in audits.values()),
        'method':'Exact rational/radical formulas and exact generator-grid parameters; saved decimals used only as an independent final check.',
        'entriesAudit':audits}
(ROOT/'audit/exact-sector-values-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='entriesAudit'},indent=2))
