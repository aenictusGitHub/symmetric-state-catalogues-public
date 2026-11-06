#!/usr/bin/env python3
"""Reproducible, actual-state witnesses for the N=3,4,5 sector atlas.

No interpolation or convex hull is used to invent a state witness.
The Pauli-normalized sectors are recovered independently from reduced purities.
"""
from pathlib import Path
from math import comb, sqrt, sin, cos, acos, asin, pi
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "data"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20261106)


def normalize(c):
    c = np.asarray(c, dtype=complex)
    c /= np.linalg.norm(c)
    pivot = int(np.argmax(np.abs(c)))
    c *= np.exp(-1j * np.angle(c[pivot]))
    c.real[np.abs(c.real) < 2e-15] = 0
    c.imag[np.abs(c.imag) < 2e-15] = 0
    return c


def sectors(c):
    c = normalize(c)
    n = len(c)-1
    s = [1.0]
    for t in range(1, n+1):
        a = np.array([[c[i+j]*sqrt(comb(t,i)*comb(n-t,j)/comb(n,i+j))
                      for j in range(n-t+1)] for i in range(t+1)])
        rho = a @ a.conj().T
        purity = float(np.sum(np.abs(rho)**2))
        s.append(comb(n,t)*(2**t*purity-sum(comb(t,k)*s[k]/comb(n,k) for k in range(t))))
    return [0.0 if abs(v)<2e-11 else float(v) for v in s[1:]]


def vector(n, terms):
    c = np.zeros(n+1, dtype=complex)
    for k, v in terms.items(): c[k] = v
    return normalize(c)


def row(c, identifier=None, label=None, **extra):
    c = normalize(c)
    r = {"coefficients": [[round(float(v.real), 13),round(float(v.imag),13)] for v in c],
         "sectors": [round(v,12) for v in sectors(c)], **extra}
    if identifier is not None: r['id'] = identifier
    if label is not None: r['label'] = label
    return r


def ghz(n, phi):
    return vector(n, {0:cos(phi/2), n:sin(phi/2)})


def n4_witness(x, y):
    """Exact inversion of PRA113 Eq.(47)-(50) for the entire triangle."""
    Y = 3-y/2
    if abs(Y)<1e-12:
        u = 1
        v = sqrt(max(0,1-x/4))
    else:
        d = x+2*Y
        u = (x+Y*sqrt(max(0,4-d)))/d
        v = 1-x*(1-u)/(Y*(1+u))
    theta = acos(min(1,max(0,u)))
    phi = asin(min(1,max(-1,v)))
    return vector(4, {0:cos(theta/2)*cos(phi/2),
                      2:1j*sin(theta/2),4:cos(theta/2)*sin(phi/2)})


def make_curve(identifier,label,points,formula,status='attained', **kwargs):
    return dict(id=identifier,label=label,status=status,formula=formula,points=points,**kwargs)


def random_rows(n,count):
    result=[]
    for j in range(count):
        c=RNG.normal(size=n+1)+1j*RNG.normal(size=n+1)
        result.append(row(c,identifier=f'n{n}-haar-{j}',label='Pure symmetric sample',
                          status='attained',source='seeded complex Gaussian Dicke amplitudes'))
    return result


def common(n):
    return dict(n=n,dimension=1 if n==3 else 2,axes=['S1'] if n==3 else ['S1','S2'],
                normalization='S_k = sum of squared expectations of all weight-k Pauli strings; S_0=1.',
                coefficientConvention='Dicke basis |D_N^(k)> with k excitations; J_z=N/2-k.',
                provenance=[{'title':'Multiqubit monogamy relations beyond shadow inequalities',
                             'doi':'10.1103/PRA.113.012415',
                             'localSource':'PRA113_012415(2026).pdf'}])


data={}
d=common(3)
d.update(identities=['S2 = 3','S3 = 4 - S1'],range='0 ≤ S1 ≤ 3',
         rangeStatus='exact attained interval',
         knownStates=[row(vector(3,{0:1}), 'n3-product','|0⟩⊗³',latex=r'|0\rangle^{\otimes3}',
                          exactSectors=['3','3','1'],status='attained'),
                      row(ghz(3,pi/2),'n3-ghz','|GHZ(3)⟩',latex=r'|\mathrm{GHZ}(3)\rangle',
                          exactSectors=['0','3','4'],status='attained'),
                      row(vector(3,{1:1}),'n3-w','|W(3)⟩',latex=r'|W(3)\rangle',
                          exactSectors=['1/3','3','11/3'],status='attained')])
d['curves']=[make_curve('n3-entire-interval','Generalized GHZ family',
                       [row(ghz(3,phi),f'n3-ghz-{j}','|GHZ(3, φ)⟩',parameter={'phi':float(phi)})
                        for j,phi in enumerate(np.linspace(0,pi/2,201))],
                       'cos(φ/2)|D₃⁽⁰⁾⟩ + sin(φ/2)|D₃⁽³⁾⟩; S1=3 cos²φ',
                       covers='entire exact set')]
d['cloud']=random_rows(3,100)
data['3']=d

d=common(4)
d.update(identities=['S3 = 12 + S1 - 2 S2','S4 = 3 - 2 S1 + S2'],
         range='0 ≤ S1 ≤ S2 − 2 ≤ 4',rangeStatus='exact attained triangle',
         polygon=[[0,2],[4,6],[0,6]],
         knownStates=[row(vector(4,{0:1}),'n4-product','|0⟩⊗⁴',latex=r'|0\rangle^{\otimes4}',status='attained',exactSectors=['4','6','4','1']),
                      row(ghz(4,pi/2),'n4-ghz','|GHZ(4)⟩',latex=r'|\mathrm{GHZ}(4)\rangle',status='attained',exactSectors=['0','6','0','9']),
                      row(vector(4,{0:1,3:sqrt(2)}),'n4-tetra','|tetra⟩',latex=r'|\mathrm{tetra}\rangle',status='attained',exactSectors=['0','2','8','5']),
                      row(vector(4,{1:1}),'n4-w','|W(4)⟩',latex=r'|W(4)\rangle',status='attained',exactSectors=['1','3','7','4']),
                      row(vector(4,{2:1}),'n4-dicke2','|D₄⁽²⁾⟩',latex=r'|D_4^{(2)}\rangle',status='attained',exactSectors=['0','6','0','9'])])
d['curves']=[make_curve('n4-lower','Lower edge',
                        [row(n4_witness(float(x),float(x+2)),f'n4-lower-{j}','Four-qubit boundary family',parameter={'S1':float(x)})
                         for j,x in enumerate(np.linspace(0,4,151))],
                        'S2=S1+2; exact Eq.(47) inverse'),
             make_curve('n4-left','Zero-dipole edge',
                        [row(n4_witness(0,float(y)),f'n4-left-{j}','Four-qubit zero-dipole family',parameter={'S2':float(y)})
                         for j,y in enumerate(np.linspace(2,6,121))],
                        'S1=0; exact Eq.(47) inverse'),
             make_curve('n4-top','Generalized GHZ edge',
                        [row(ghz(4,float(phi)),f'n4-top-{j}','|GHZ(4, φ)⟩',parameter={'phi':float(phi)})
                         for j,phi in enumerate(np.linspace(0,pi/2,151))],
                        'S2=6; cos(φ/2)|D₄⁽⁰⁾⟩+sin(φ/2)|D₄⁽⁴⁾⟩')]
d['cloud']=[]
for j in range(51):
    x=4*j/50
    for k in range(41):
        y=x+2+(4-x)*k/40
        w=row(n4_witness(x,y),f'n4-grid-{j}-{k}','Exact four-qubit witness',status='attained',source='PRA113 Eq.(47)-(50), exact inverse')
        assert max(abs(np.array(w['sectors'][:2])-[x,y]))<2e-9
        d['cloud'].append(w)
d['cloud']+=random_rows(4,200)
d['provenance'][0]['equations']='(44)-(50)'
data['4']=d

rt5=sqrt(5)
pc=2-3/rt5
p0=rt5/3
xc=(22805-9840*rt5)/11881


def rho(p):
    return max(0,2*p*(5*p*p+6*rt5*p-11)/(3*rt5*p*p+2*p+3*rt5))


def n5_small(p):
    return vector(5,{0:-p,2:sqrt(rho(p)),4:1})


def n5_root(lam):
    return vector(5,{1:1,5:lam})


def n5_left(t):
    return vector(5,{0:1,2:sqrt(max(0,3*t*t-5)),4:t})


d=common(5)
d.update(identities=['S3 = 10 + 2 S1 - S2','S4 = 15 - S2','S5 = 6 - 3 S1 + S2'],
         range='0 ≤ x ≤ 5, b(x) ≤ y ≤ 10',rangeStatus='exact range in dedicated completed five-qubit proof',
         rigorousLinearOuter='0 ≤ S1 ≤ 5; 10/3 ≤ S2 ≤ 10; 2 S1 ≤ S2',
         lowerBoundary={'junctionX':xc,'junctionXExact':'(22805−9840√5)/11881',
                        'pieces':['Parametric algebraic C2 arc for 0≤x≤xc',
                                  '(3x−2√(5x)+15)/4 for xc≤x≤9/5',
                                  '2x for 9/5≤x≤5']},
         proofStatusNote='Uses the dedicated completed proof, not the older arxiv-161026.tex text where this reduction is still marked open.',
         knownStates=[row(vector(5,{0:1}),'n5-product','|0⟩⊗⁵',latex=r'|0\rangle^{\otimes5}',status='attained',exactSectors=['5','10','10','5','1']),
                      row(ghz(5,pi/2),'n5-ghz','|GHZ(5)⟩',latex=r'|\mathrm{GHZ}(5)\rangle',status='attained',exactSectors=['0','10','0','5','16']),
                      row(vector(5,{1:1}),'n5-w','|W(5)⟩',latex=r'|W(5)\rangle',status='attained',exactSectors=['9/5','18/5','10','57/5','21/5']),
                      row(vector(5,{2:1}),'n5-dicke2','|D₅⁽²⁾⟩',latex=r'|D_5^{(2)}\rangle',status='attained'),
                      row(vector(5,{1:sqrt(5),5:1}),'n5-psia','|ψ₅⁽ᵃ⁾⟩',latex=r'|\psi_5^{(a)}\rangle',status='attained',exactSectors=['5/9','10/3','70/9','35/3','23/3']),
                      row(vector(5,{0:rt5,2:sqrt(2),4:-3}),'n5-psib','|ψ₅⁽ᵇ⁾⟩',latex=r'|\psi_5^{(b)}\rangle',status='attained',exactSectors=['0','18/5','32/5','57/5','48/5']),
                      row(n5_root(pc),'n5-psic','|ψ₅⁽ᶜ⁾⟩',latex=r'|\psi_5^{(c)}\rangle',status='attained')])
d['curves']=[make_curve('n5-small-arc','Algebraic C₂ arc',
                       [row(n5_small(float(p)),f'n5-small-{j}','C₂ boundary state',parameter={'p':float(p)},symmetry='C2')
                        for j,p in enumerate(np.linspace(p0,pc,121))],
                       '(-p D₀ + √ρ(p) D₂ + D₄)/√(1+p²+ρ(p))'),
             make_curve('n5-root-arc','Square-root C₄ arc',
                       [row(n5_root(float(lam)),f'n5-root-{j}','C₄ boundary state',parameter={'lambda':float(lam)},symmetry='C4')
                        for j,lam in enumerate(np.linspace(pc,0,181))],
                       '(D₁+λ D₅)/√(1+λ²)'),
             make_curve('n5-linear-arc','Product–W boundary',
                       [row(vector(5,{0:cos(theta),1:sin(theta)}),f'n5-line-{j}','Product–W boundary state',parameter={'theta':float(theta)})
                        for j,theta in enumerate(np.linspace(pi/2,0,181))],
                       'cosθ D₀+sinθ D₁; S2=2S1'),
             make_curve('n5-top','Generalized GHZ edge',
                       [row(ghz(5,float(phi)),f'n5-top-{j}','|GHZ(5, φ)⟩',parameter={'phi':float(phi)})
                        for j,phi in enumerate(np.linspace(0,pi/2,181))],
                       'cos(φ/2)D₀+sin(φ/2)D₅; S2=10'),
             make_curve('n5-left-negative','Zero-dipole family, lower part',
                       [row(n5_left(float(t)),f'n5-left-neg-{j}','C₂ zero-dipole state',parameter={'t':float(t)})
                        for j,t in enumerate(np.linspace(-3/rt5,-sqrt((345+60*sqrt(15))/289),101))],
                       '(D₀+√(3t²−5)D₂+tD₄)/(2√(t²−1)); S1=0'),
             make_curve('n5-left-positive','Zero-dipole family, upper part',
                       [row(n5_left(float(t)),f'n5-left-pos-{j}','C₂ zero-dipole state',parameter={'t':float(t)})
                        for j,t in enumerate(np.linspace(sqrt(5/3),rt5,181))],
                       '(D₀+√(3t²−5)D₂+tD₄)/(2√(t²−1)); S1=0')]
d['cloud']=random_rows(5,2400)
for j in range(1300):
    # Parity and sparse families resolve the small-S1 portion much more clearly
    # than Haar-only samples. Every point still carries its actual state.
    c=np.zeros(6,dtype=complex)
    support = [0,2,4] if j%2 else RNG.choice(6,size=2,replace=False)
    c[support]=RNG.normal(size=len(support))+1j*RNG.normal(size=len(support))
    d['cloud'].append(row(c,f'n5-parity-{j}','Sparse symmetric sample',status='attained'))
d['provenance'] += [
    {'title':'The exact (S1,S2) range of pure symmetric five-qubit states: A compact rigorous proof',
     'localSource':'output/pdf/N5_symmetric_S1S2_boundary_compact_proof.tex',
     'section':'Theorem 1.2, Eq.(main-set), state parametrizations (state-v2), (state-infty), (state-line)'},
    {'title':'Dedicated manuscript with completed N5 symmetry reduction',
     'localSource':'arxiv-161026proofN5.tex',
     'section':'Theorems N5C2scope and N5exactLowerBoundary'}]
d['polygon'] = [p['sectors'][:2] for curve in d['curves'][:4] for p in curve['points']]
d['polygon'] += [p['sectors'][:2] for curve in reversed(d['curves'][4:]) for p in reversed(curve['points'])]
d['polygonNote'] = 'Boundary tessellation of the exact proven two-dimensional range; all displayed boundary samples have state witnesses. No convex hull of samples is used.'
data['5']=d

# Validate every output state and all pure-state identities before serializing.
maxerr=0.0
counts={}
for key,d in data.items():
    n=int(key)
    records=d['knownStates']+d['cloud']+[p for curve in d['curves'] for p in curve['points']]
    for r in records:
        c=np.array([complex(*v) for v in r['coefficients']])
        assert abs(np.vdot(c,c).real-1)<2e-12
        s=np.array(sectors(c))
        maxerr=max(maxerr,float(np.max(np.abs(s-np.array(r['sectors'])))))
        if 'exactSectors' in r:
            from fractions import Fraction
            expected=np.array([float(Fraction(v)) for v in r['exactSectors']])
            assert max(abs(s-expected))<1e-9,(n,r['id'],s,expected)
        assert abs(1+sum(s)-2**n)<1e-9
        if n==3: assert max(abs(s-np.array([s[0],3,4-s[0]])))<1e-9
        if n==4: assert max(abs(s[2:]-np.array([12+s[0]-2*s[1],3-2*s[0]+s[1]])))<1e-9
        if n==5: assert max(abs(s[2:]-np.array([10+2*s[0]-s[1],15-s[1],6-3*s[0]+s[1]])))<1e-9
    counts[key]={'states':len(d['knownStates']),'curves':len(d['curves']),'cloud':len(d['cloud']),'witnesses':len(records)}

payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
(OUT/'small-n.json').write_text(payload+'\n')
(OUT/'small-n.js').write_text('/* Generated actual-state witnesses; see generate_small_n.py. */\n'+'globalThis.SECTOR_SMALL_DATA = '+payload+';\n')
audit={'status':'PASS','seed':20261106,'counts':counts,'maxSectorRoundtripError':maxerr,
       'method':'Dicke bipartition amplitudes → reduced purities → Pauli sector binomial inversion',
       'n4GridInversionChecked':True,'allNamedRationalSectorsChecked':True}
(OUT/'small-n-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
