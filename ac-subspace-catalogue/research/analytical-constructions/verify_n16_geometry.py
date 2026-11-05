#!/usr/bin/env python3
"""Fast exact multiplicity certificates and coordinates for the N16 t2 basis.

Avoids constructing a primitive element for a multiquadratic coefficient field.
The reciprocal quartic is square-free iff its two values at +/-1 and its
quadratic discriminant in y+1/y are nonzero. The other reduced polynomials
are cubics, checked by their explicit discriminants. Nonzero imaginary parts
are certified with rational enclosing bounds on integer square roots.

Use --write to generate n16-geometry.json. Without it, verify the saved
coordinates/certificates against the current exact coefficients. Requires
SymPy, NumPy and mpmath. All paths are relative to this script.
"""
from pathlib import Path
from math import comb,isqrt,sqrt
import argparse,cmath,hashlib,json,math
import sympy as s
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parent
SPACE='mod4-reflection-n16-d3-t2'
y=s.Symbol('y')

def simple(x):return s.simplify(s.expand_complex(x))
def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()

def interval(expression):
    expression=s.expand(expression);lo=s.S.Zero;hi=s.S.Zero;terms=[]
    for term in s.Add.make_args(expression):
        coefficient,rest=term.as_coeff_Mul()
        assert coefficient.is_Rational
        if rest==1:lower=upper=s.S.One;radicand=None
        else:
            assert rest.is_Pow and rest.exp==s.Rational(1,2) and rest.base.is_Integer and rest.base>0,rest
            radicand=int(rest.base);floor=isqrt(radicand);lower=s.Integer(floor)
            upper=lower if floor*floor==radicand else lower+1
        assert radicand is None or lower*lower<=radicand<=upper*upper
        lo+=coefficient*(lower if coefficient>=0 else upper)
        hi+=coefficient*(upper if coefficient>=0 else lower)
        terms.append({'coefficient':str(coefficient),'radicand':radicand,'sqrtLower':str(lower),'sqrtUpper':str(upper)})
    return lo,hi,terms

def nonzero(name,expression):
    expression=simple(expression)
    for part,component in [('imaginary',simple(s.im(expression))),('real',simple(s.re(expression)))]:
        if component==0:continue
        lo,hi,terms=interval(component)
        if hi<0 or lo>0:
            return {'name':name,'expression':str(expression),'nonzeroPart':part,
                    'component':str(component),'rationalInterval':[str(lo),str(hi)],
                    'intervalTerms':terms,'nonzeroVerifiedExactly':True}
    raise AssertionError(('No nonzero interval established',name,expression))

def polynomial_certificate(raw,index):
    coefficients=list(map(s.sympify,raw['exactCoefficients']));support=[r for r,c in enumerate(coefficients)if c!=0]
    lo,hi=min(support),max(support)
    assert support==list(range(lo,hi+1,4))
    terms=[simple((-1)**r*s.sqrt(s.binomial(16,r))*coefficients[r])for r in support]
    leading=terms[-1];terms=[simple(c/leading)for c in terms]
    checks=[nonzero('constant coefficient',terms[0])]
    if index==0:
        assert len(terms)==5 and terms[0]==terms[4]==1 and terms[1]==terms[3]
        a,b=terms[1:3]
        checks.extend([nonzero('q(1)=b+2+2a',b+2+2*a),nonzero('q(-1)=b+2-2a',b+2-2*a),nonzero('quadratic discriminant a^2-4(b-2)',a*a-4*(b-2))])
        method='Reciprocal quartic: q(y)=y^2[(y+1/y)^2+a(y+1/y)+(b-2)]; nonzero q(+/-1) and distinct quadratic roots prove square-freeness.'
    else:
        assert len(terms)==4 and terms[-1]==1
        d,c,b,a=terms
        disc=b*b*c*c-4*a*c*c*c-4*b*b*b*d-27*a*a*d*d+18*a*b*c*d
        checks.append(nonzero('cubic discriminant',disc))
        method='Explicit cubic discriminant is nonzero, so all three nonzero roots are simple.'
    certificate={'method':method,'reducedPolynomialAscending':list(map(str,terms)),
       'nonzeroChecks':checks,'zeroMultiplicity':lo,'infinityMultiplicity':16-hi,'substitutionStep':4,
       'conclusion':'All nonzero finite roots of P(z)=z^low q(z^4) are simple; only the explicit zero and infinity factors can repeat.'}
    return coefficients,terms,certificate

def mpc(value):
    return mp.mpc(str(s.N(s.re(value),90)),str(s.N(s.im(value),90)))

def point(z,m=1):
    if z is None:return {'x':0.,'y':0.,'z':1.,'multiplicity':m}
    q=1+abs(z)**2
    return {'x':float(2*z.real/q),'y':float(-2*z.imag/q),'z':float((abs(z)**2-1)/q),'multiplicity':m}

def roots_and_points(terms,certificate):
    mp.mp.dps=80;values=list(map(mpc,terms))
    if len(values)==5:
        a,b=values[1:3];delta=mp.sqrt(a*a-4*(b-2));roots=[]
        for t in [(-a+delta)/2,(-a-delta)/2]:
            branch=mp.sqrt(t*t-4);roots.extend([(t+branch)/2,(t-branch)/2])
    else:roots=mp.polyroots(list(reversed(values)),maxsteps=500,error=False)
    residual=max(abs(sum(c*r**i for i,c in enumerate(values)))/sum(abs(c*r**i)for i,c in enumerate(values))for r in roots)
    assert residual<mp.mpf('1e-65')
    points=[]
    if certificate['zeroMultiplicity']:points.append(point(mp.mpc(0),certificate['zeroMultiplicity']))
    if certificate['infinityMultiplicity']:points.append(point(None,certificate['infinityMultiplicity']))
    for value in roots:
        radius=abs(value)**mp.mpf('.25');angle=mp.arg(value)
        for k in range(4):points.append(point(mp.rect(radius,(angle+2*mp.pi*k)/4)))
    return points,mp.nstr(residual,20)

def ray_error(points,exact):
    poly=np.array([1.+0j])
    for p in points:
        a=sqrt(max(0.,(1+p['z'])/2))
        b=(p['x']+1j*p['y'])/(2*a)if a>1e-10 else 1.+0j
        for _ in range(p['multiplicity']):poly=np.convolve(poly,[a,-b])
    restored=np.array([(-1)**r*c/sqrt(comb(16,r))for r,c in enumerate(poly)])
    expected=np.array([complex(s.N(c,30))for c in exact])
    expected/=np.linalg.norm(expected);restored/=np.linalg.norm(restored)
    overlap=np.vdot(expected,restored)
    return float(np.linalg.norm(expected-restored*np.conj(overlap)/abs(overlap)))

def check_points(points,exact,certificate):
    assert sum(p['multiplicity']for p in points)==16
    expected=sorted([certificate['zeroMultiplicity'],certificate['infinityMultiplicity']]+[1]*(4*(len(certificate['reducedPolynomialAscending'])-1)),reverse=True)
    expected=[x for x in expected if x]
    assert sorted([p['multiplicity']for p in points],reverse=True)==expected
    unit=max(abs(sum(p[key]**2 for key in ['x','y','z'])-1)for p in points)
    separation=min(math.sqrt(sum((a[key]-b[key])**2 for key in ['x','y','z']))for i,a in enumerate(points)for b in points[i+1:])
    error=ray_error(points,exact)
    assert unit<1e-12 and separation>1e-7 and error<1e-10
    return {'projectiveDickeReconstructionError':error,'maxUnitSphereError':unit,'minimumDistinctPointDistance':separation},expected

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    source=json.loads((ROOT/'t2-bases.json').read_text());space=next(st for st in source['subspaces']if st['id']==SPACE)
    output_path=ROOT/'n16-geometry.json'
    previous=None if args.write else json.loads(output_path.read_text())
    records={};audits=[]
    for i,raw in enumerate(space['basis']):
        bid=SPACE+'-basis-'+str(i+1);exact,terms,certificate=polynomial_certificate(raw,i)
        if args.write:points,residual=roots_and_points(terms,certificate)
        else:
            stored=previous['records'][bid]
            assert stored['sourceExactCoefficients']==raw['exactCoefficients'] and stored['sourceCoefficientsSHA256']==digest(raw['exactCoefficients'])
            assert stored['diagnostics']['exactSquareFreeCertificate']==certificate
            points=stored['points'];residual=stored['diagnostics']['maximumRelativePolynomialResidual']
        checked,partition=check_points(points,exact,certificate)
        saved=np.array([complex(*z) for z in raw['coefficients']]);exact_float=np.array([complex(s.N(c,30))for c in exact])
        saved/=np.linalg.norm(saved);exact_float/=np.linalg.norm(exact_float)
        phase=np.vdot(exact_float,saved)
        formula_error=float(np.linalg.norm(exact_float-saved*np.conj(phase)/abs(phase)))
        assert formula_error<1e-12
        records[bid]={'id':bid,'n':16,'points':points,'partition':partition,
            'sourceExactCoefficients':raw['exactCoefficients'],'sourceCoefficientsSHA256':digest(raw['exactCoefficients']),
            'coordinateMethod':'Exact elementary square-freeness certificates and 80-digit root placement; no primitive-element algebra or numerical clustering.',
            'multiplicityEvidence':'Exact nonvanishing certificates for the quartic reduction or cubic discriminant.',
            'diagnostics':{'coefficientSource':'saved exact Dicke coefficients','squareFreeFactorizationExact':True,
                'zeroMultiplicity':certificate['zeroMultiplicity'],'infinityMultiplicity':certificate['infinityMultiplicity'],'nonzeroRootStep':4,
                'reducedPolynomialFactors':[{'polynomial':str(sum(c*y**r for r,c in enumerate(terms))),'multiplicity':1}],
                'exactSquareFreeCertificate':certificate,'rootPrecisionDecimalDigits':80,'maximumRelativePolynomialResidual':residual,
                'exactFormulaVsSavedCoefficientsError':formula_error,**checked}}
        audits.append({'id':bid,'exactMultiplicityCertificateVerified':True,'partition':partition,**checked})
        print('Verified',bid,'ray error',checked['projectiveDickeReconstructionError'],flush=True)
    output={'schemaVersion':1,'source':'t2-bases.json','subspaceId':SPACE,'records':records,
            'convention':'P(z)=sum_r (-1)^r sqrt(C(16,r)) c_r z^r; root z maps to (2Re z,-2Im z,|z|^2-1)/(1+|z|^2), infinity north.',
            'allMultiplicityCertificatesVerifiedExactly':True,'noNumericalRootClustering':True}
    if args.write:output_path.write_text(json.dumps(output,indent=2)+'\n')
    report={'allChecksPassed':True,'method':'Exact reciprocal-quartic/cubic-discriminant tests, rational square-root bounds and independent Dicke reconstruction.',
            'records':audits,'maximumProjectiveDickeReconstructionError':max(x['projectiveDickeReconstructionError']for x in audits)}
    (ROOT/'n16-geometry-verification.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
