#!/usr/bin/env python3
"""Build offline Majorana coordinates from this catalogue's Dicke representatives.

Requires Python 3, SymPy and NumPy. Run from any directory. The original
catalogue is read only; this creates data/majorana.{json,js}.

P(z) = sum_r (-1)^r sqrt(binomial(n,r)) c_r z^r.
A finite root z corresponds to the spinor proportional to (z,1), hence
Bloch point (2 Re(z), -2 Im(z), |z|^2-1)/(1+|z|^2).
The projective root at infinity is the north pole. Repeated roots are
combined only when exact square-free algebra establishes multiplicity;
no distance-based clustering is used.
"""
from pathlib import Path
from functools import reduce
from math import gcd, comb, pi
import cmath
import hashlib
import json
import math
import re
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
y = s.Symbol('y')

# Exact formula representatives displayed in the sector catalogue. These are
# independently compared with every saved Dicke vector before they are used.
SECTOR_FORMULAS = {
 'tetra': {0:'1/sqrt(3)',3:'sqrt(2)/sqrt(3)'},
 'psi5a': {1:'sqrt(5/6)',5:'sqrt(1/6)'},
 'octa': {1:'1/sqrt(2)',5:'1/sqrt(2)'},
 'psi6': {2:'sqrt(20)/sqrt(21)',5:'1/sqrt(21)'},
 'pyramid': {0:'-1/2',4:'sqrt(3)/2'},
 'psi7': {0:'2*sqrt(2)/6',3:'sqrt(14)/6',6:'sqrt(14)/6'},
 'psi7b': {1:'-1/sqrt(18)',4:'sqrt(10)/sqrt(18)',7:'sqrt(7)/sqrt(18)'},
 'psi8': {0:'sqrt(30)/12',4:'2*sqrt(21)/12',8:'sqrt(30)/12'},
 'psi8b': {2:'4/(3*sqrt(3))',5:'2/(3*sqrt(3))',8:'sqrt(7)/(3*sqrt(3))'},
 'tetra8': {0:'sqrt(7)/(3*sqrt(3))',3:'2/(3*sqrt(3))',6:'4/(3*sqrt(3))'},
 'psi9': {0:'1/sqrt(3)',3:'sqrt(2)/sqrt(3)'},
 'psi9b': {0:'5/sqrt(91)',6:'sqrt(66)/sqrt(91)'},
 'psi10': {1:'1/sqrt(2)',7:'1/sqrt(2)'},
 'psi10c': {1:'sqrt(5)/4',5:'sqrt(6)/4',9:'sqrt(5)/4'},
 'psi10d': {0:'1/2',5:'sqrt(2)/2',10:'1/2'},
 'psi11': {2:'sqrt(17)/sqrt(18)',8:'1/sqrt(18)'},
 'psi12': {2:'sqrt(7)/sqrt(8)',6:'1/sqrt(8)'},
}

def parse(text):
    # rational=True prevents a displayed rational such as 5/6 becoming a float.
    return s.sympify(str(text), rational=True)

def as_complex(value):
    return complex(*value) if isinstance(value,list) else complex(value)

def normed(values):
    values=np.array(values,dtype=complex)
    return values/np.linalg.norm(values)

def ray_error(a,b):
    a=normed(a); b=normed(b)
    overlap=np.vdot(a,b)
    return float(np.linalg.norm(a-b*np.conj(overlap)/abs(overlap)))

def exact_coefficients(record):
    n=record['n']; source=record.get('exactCoefficients')
    if source:
        if isinstance(source,dict):
            return [parse(source.get(str(r),'0')) for r in range(n+1)],'saved exact Dicke coefficients'
        return list(map(parse,source)),'saved exact Dicke coefficients'
    key=record.get('stateKey')
    values=None
    if key=='coherent': values={0:1}
    elif key=='ghz':values={0:'1/sqrt(2)',n:'1/sqrt(2)'}
    elif key=='dicke':
        support=[r for r,c in enumerate(record['coefficients']) if abs(as_complex(c))>0]
        if len(support)==1:values={support[0]:1}
    elif key=='psi9e':
        t=(83+9*s.sqrt(669))/430
        u=(37052+171*s.sqrt(669))/210700
        denominator=s.sqrt(1+14*u+4*t*t)
        values={1:1/denominator,4:s.sqrt(14*u)/denominator,7:2*t/denominator}
    elif key in SECTOR_FORMULAS: values=SECTOR_FORMULAS[key]
    if values is not None:
        return [parse(values.get(r,0)) for r in range(n+1)],'exact formula displayed in the catalogue'
    return None,'saved decimal Dicke coefficients'

def point(root,m):
    if root is None:return {'x':0.,'y':0.,'z':1.,'multiplicity':int(m)}
    q=1+abs(root)**2
    return {'x':2*root.real/q,'y':-2*root.imag/q,'z':(abs(root)**2-1)/q,'multiplicity':int(m)}

def reconstruct(points,n):
    poly=np.array([1.+0j])
    for p in points:
        z=p['z']; a=math.sqrt(max(0.,(1+z)/2))
        b=(p['x']+1j*p['y'])/(2*a) if a>1e-10 else 1.+0j
        for _ in range(p['multiplicity']):poly=np.convolve(poly,[a,-b])
    return np.array([(-1)**r*poly[r]/math.sqrt(comb(n,r)) for r in range(n+1)])

CACHE={}
IMPORTED_GEOMETRY=json.loads((ROOT/'research/cross-catalogue/exact-states.json').read_text()).get('geometryImports',{})
def coordinates(record):
    n=record['n']; exact,source=exact_coefficients(record)
    saved=[as_complex(c) for c in record['coefficients']]
    if record.get('status')=='numerical' and record.get('majoranaPoints'):
        points=record['majoranaPoints']
        assert sum(p['multiplicity'] for p in points)==n
        error=ray_error(saved,reconstruct(points,n));assert error<1e-10,(record['id'],error)
        return {'id':record['id'],'n':n,'label':record['label'],'points':points,
          'partition':record['partition'],'multiplicityEvidence':'Imported certified partition of the saved decimal vector; AC remains numerical',
          'evidence':f"Numerical library representative, reported AC order {record['achievedOrder']}. Reduced-state residuals are checked from the saved 85-digit components. This is not an exact AC or optimality proof.",
          'coordinateMethod':'Imported audited refined Majorana points, independently reconstructed against the displayed Dicke vector.',
          'diagnostics':{'coefficientSource':'Saved 85-digit numerical library state; browser coefficients rounded','squareFreeFactorizationExact':False,'projectiveDickeReconstructionError':error,'exactFormulaVsSavedCoefficientsError':0.,'maxUnitSphereError':max(abs(sum(p[k]**2 for k in ['x','y','z'])-1) for p in points),'noNumericalRootClustering':True}}
    if record['id'] in IMPORTED_GEOMETRY:
        imported=IMPORTED_GEOMETRY[record['id']]
        assert imported['exactCoefficients']==record['exactCoefficients']
        points=imported['points'];diagnostics=dict(imported['diagnostics'])
        assert diagnostics['squareFreeFactorizationExact'] and sum(p['multiplicity'] for p in points)==n
        error=ray_error(saved,reconstruct(points,n));assert error<1e-10
        diagnostics['projectiveDickeReconstructionError']=error
        diagnostics['importedExactGeometryCertificate']=imported['sourceGeometryCertificate']
        return {'id':record['id'],'n':n,'label':record['label'],'points':points,'diagnostics':diagnostics,
          'evidence':f"Exact witness of AC order {record['achievedOrder']}. Coordinates render its certified radical Dicke coefficients.",
          'coordinateMethod':'Imported exact square-free factorization and numerical roots; the state is independently reconstructed from these points.',
          'partition':sorted([p['multiplicity'] for p in points],reverse=True),'multiplicityEvidence':'exact polynomial factorization'}
    if exact is not None:
        exact_numeric=[complex(s.N(c,40)) for c in exact]
        formula_error=ray_error(saved,exact_numeric)
        assert formula_error<3e-12,(record['id'],'exact formula differs',formula_error)
        coefficients=exact
    else:
        coefficients=[s.Rational(str(c.real))+s.I*s.Rational(str(c.imag)) for c in saved]
        formula_error=0.
    key=(n,tuple(map(str,coefficients)),exact is not None)
    if key not in CACHE:
        support=[r for r,c in enumerate(coefficients) if c!=0]
        low,high=min(support),max(support)
        step=reduce(gcd,(r-low for r in support),0) or 1
        values={int((r-low)//step):s.simplify((-1)**r*s.sqrt(s.binomial(n,r))*coefficients[r]) for r in support}
        common=values[0]
        # Removing a common scalar keeps the algebraic number field small.
        poly=s.Poly(sum(s.simplify(c/common)*y**r for r,c in values.items()),y,extension=True)
        if exact is not None:
            _,factors=s.sqf_list(poly)
        else:
            # Numerical states are never clustered or asserted to have exact
            # finite repeated roots. The explicit zero/infinity support stays.
            factors=[(poly,1)]
        points=[]
        if low:points.append(point(0j,low))
        if high<n:points.append(point(None,n-high))
        descriptions=[]
        for factor,m in factors:
            descriptions.append({'polynomial':str(factor.as_expr()),'multiplicity':int(m)})
            roots=s.nroots(factor,n=40,maxsteps=500)
            for root in roots:
                root=complex(root); radius=abs(root)**(1/step); angle=cmath.phase(root)
                for k in range(step):points.append(point(cmath.rect(radius,(angle+2*pi*k)/step),m))
        assert sum(p['multiplicity'] for p in points)==n
        restored=reconstruct(points,n)
        reconstruction_error=ray_error(saved,restored)
        assert reconstruction_error<1e-10,(record['id'],reconstruction_error)
        min_separation=min((math.sqrt(sum((a[k]-b[k])**2 for k in ['x','y','z'])) for i,a in enumerate(points) for b in points[i+1:]),default=None)
        if min_separation is not None: assert min_separation>1e-7,(record['id'],'unresolved close roots',min_separation)
        max_unit_error=max(abs(sum(p[k]**2 for k in ['x','y','z'])-1) for p in points)
        CACHE[key]={'points':points,'diagnostics':{
            'coefficientSource':source,'squareFreeFactorizationExact':exact is not None,
            'zeroMultiplicity':low,'infinityMultiplicity':n-high,'nonzeroRootStep':step,
            'reducedPolynomialFactors':descriptions,
            'projectiveDickeReconstructionError':reconstruction_error,
            'exactFormulaVsSavedCoefficientsError':formula_error,
            'maxUnitSphereError':max_unit_error,'minimumDistinctPointDistance':min_separation,
            'rootPrecisionDecimalDigits':40,
        }}
    result=json.loads(json.dumps(CACHE[key]))
    label=record.get('label') or record.get('stateName') or record.get('state') or 'Representative state'
    if 'achievedOrder' in record:
        evidence=f"Exact witness of AC order {record['achievedOrder']}. Coordinates are a numerical rendering of its exact Dicke coefficients."
    elif record.get('status')=='numerical':
        evidence='Numerical candidate. The sphere represents the saved Dicke coefficients; global optimality is not proved.'
    elif record.get('status')=='proven':
        evidence='Representative of the proved sector-length extremum. Coordinates are a numerical rendering of the saved state.'
    else:
        evidence='Representative of a candidate sector-length extremum. Global optimality is not proved for this entry.'
    if exact is None:
        evidence+=' Coordinates use the saved decimal coefficients; coincident finite roots are not inferred by clustering.'
    result.update(id=record['id'],n=n,label=label,evidence=evidence,
       coordinateMethod='Exact square-free factorization followed by numerical roots.' if exact is not None else 'Numerical roots of the saved decimal polynomial, retaining its explicit zero and infinity factors.',
       partition=sorted([p['multiplicity'] for p in result['points']],reverse=True),
       multiplicityEvidence='exact polynomial factorization' if exact is not None else 'numerical representation only')
    return result

def main():
    raw=(ROOT/'data/catalogue.json').read_bytes(); data=json.loads(raw)
    inputs=data.get('states',data['entries'])
    records={}
    for st in inputs:
        print(st['id'],flush=True)
        records[st['id']]=coordinates(st)
    output={'schemaVersion':1,'catalogueSHA256':hashlib.sha256(raw).hexdigest(),
        'convention':{'polynomial':'P(z)=sum_r (-1)^r sqrt(binomial(n,r)) c_r z^r',
                      'sphere':'(2 Re(z), -2 Im(z), |z|^2-1)/(1+|z|^2); infinity is north',
                      'basis':'Normalized Dicke basis, r counts excitations.',
                      'rotation':'Coordinates follow the saved representative; viewer rotations change only the view.'},
        'records':records,
        'validation':{'recordCount':len(records),'allRecordsHaveCoordinates':len(records)==len(inputs),
            'exactMultiplicityRecords':sum(r['diagnostics']['squareFreeFactorizationExact'] for r in records.values()),
            'maximumProjectiveDickeReconstructionError':max(r['diagnostics']['projectiveDickeReconstructionError'] for r in records.values()),
            'maximumUnitSphereError':max(r['diagnostics']['maxUnitSphereError'] for r in records.values()),
            'noNumericalRootClustering':True}}
    (ROOT/'data/majorana.json').write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n')
    (ROOT/'data/majorana.js').write_text('window.MAJORANA_DATA = '+json.dumps(output,separators=(',',':'),ensure_ascii=False)+';\n')
    index=ROOT/'index.html'
    index.write_text(re.sub(r'src="data/majorana\.js(?:\?v=[^"]*)?"','src="data/majorana.js?v='+hashlib.sha256((ROOT/'data/majorana.js').read_bytes()).hexdigest()[:12]+'"',index.read_text()))
    print(json.dumps(output['validation'],indent=2))

if __name__=='__main__':main()
